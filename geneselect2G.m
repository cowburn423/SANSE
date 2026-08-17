function [tfg, C] = geneselect2G(Ngen, Gmax, C, yERR, yOBS, CALC, maskCALC)
    % gemini modiified for gpu rand calls  Aug 26
    % Execute Ngen generations for Gmax genes natively on a single GPU
    u_mask = (nargin >= 7 && ~isempty(maskCALC));
    usegpu = canUseGPU();

    if u_mask
        if numel(maskCALC) ~= size(CALC, 1)
            error('geneselect2: bad mask size');
        end
        XR = int32(find(~maskCALC));
        if usegpu
            % Shuffle XR once on GPU via random sort
            [~, rand_idx] = sort(gpuArray.rand(numel(XR), 1, 'single'));
            XR_dev = gpuArray(XR(rand_idx));
        else
            XR_dev = XR(randperm(numel(XR)));
        end
        n_pool = numel(XR_dev);
    else
        n_pool = size(CALC, 1);
    end

    % Pre-allocate arrays in GPU device memory
    if usegpu
        CALC_dev = gpuArray(single(CALC));
        yERR_dev = gpuArray(single(yERR));
        yOBS_dev = gpuArray(single(yOBS));
        IP_mat   = gpuArray(reshape(int32([C.ip]), Gmax, numel(C)));
        oldtf    = gpuArray(single([C.tfg]));
    else
        CALC_dev = single(CALC);
        yERR_dev = single(yERR);
        yOBS_dev = single(yOBS);
        IP_mat   = reshape(int32([C.ip]), Gmax, numel(C));
        oldtf    = single([C.tfg]);
    end

    numC = numel(C);
    n_delta = floor(numC / Gmax);

    [oldtf, inC] = sort(oldtf);
    IP_mat = IP_mat(:, inC);
    OLDMIN = oldtf(1);

    tfg_history = zeros(numC, Ngen, 'single');
    
    % Pre-build coordinate grids for matrix indexing on device
    col_idx_all = 1:numC;
    row_grid = (1:Gmax)';
    if usegpu
        col_idx_all = gpuArray(col_idx_all);
        row_grid    = gpuArray(row_grid);
    end

    for igen = 1:Ngen
        [oldtf, inC] = sort(oldtf);
        IP_mat = IP_mat(:, inC);

        if oldtf(1) < OLDMIN
            OLDMIN = oldtf(1);
        else
            n_delta = n_delta + 1;
        end

        % Clone elite chromosome across all columns
        best_ip = IP_mat(:, 1);
        IP_mat(:, 2:numC) = repmat(best_ip, 1, numC - 1);

        % Vectorized stage counts (1 x numC)
        istages = min(max(1, floor(col_idx_all / n_delta)), Gmax);

        if usegpu
            % 1. GPU random replacement positions (Gmax x numC)
            [~, pos_perm] = sort(gpuArray.rand(Gmax, numC, 'single'), 1);

            % 2. GPU random candidate gene sampling
            if n_pool <= 50000
                % Exact sampling without replacement via sort
                [~, gene_perm] = sort(gpuArray.rand(n_pool, numC, 'single'), 1);
                raw_genes = int32(gene_perm(1:Gmax, :));
            else
                % Fast O(1) integer sampling on GPU (negligible collisions when Gmax << n_pool)
                raw_genes = gpuArray.randi([1, n_pool], Gmax, numC, 'int32');
            end

            if u_mask
                candidate_genes = XR_dev(raw_genes);
            else
                candidate_genes = raw_genes;
            end

            % 3. Apply mutations in-place across entire matrix in one kernel launch
            mask_mat = (row_grid <= istages);
            % Exclude elite column (col 1) from mutation
            mask_mat(:, 1) = false;

            col_grid = repmat(col_idx_all, Gmax, 1);
            mutation_indices = sub2ind([Gmax, numC], pos_perm(mask_mat), col_grid(mask_mat));
            IP_mat(mutation_indices) = candidate_genes(mask_mat);
        else
            % CPU fallback
            for ii = 2:numC
                istage = istages(ii);
                ipos = randperm(Gmax, istage);
                if u_mask
                    stage = XR_dev(randperm(n_pool, istage));
                else
                    stage = int32(randperm(n_pool, istage));
                end
                IP_mat(ipos, ii) = stage;
            end
        end

        % Unpack into struct array for target function evaluation
        IP_host = gather(IP_mat);
        for jj = 1:numC
            C(jj).ip = IP_host(:, jj);
        end

        % Target function evaluation on GPU
        current_tfg = gettf(yERR_dev, C, yOBS_dev, CALC_dev);
        if isa(current_tfg, 'gpuArray')
            oldtf = current_tfg(:)';
        else
            oldtf = gpuArray(single(current_tfg(:)'));
        end

        tfg_history(:, igen) = gather(oldtf);

        oldtf_host = gather(oldtf);
        for jj = 1:numC
            C(jj).tfg = oldtf_host(jj);
        end
    end

    [~, inC] = sort(oldtf_host);
    C = C(inC);
    tfg = tfg_history;
end