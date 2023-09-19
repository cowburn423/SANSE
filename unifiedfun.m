function Cfin= pgenslect(Ngen , Gmax, C2, yERR, yOBS, CALC, maskCALC)
% parallelise genselect2 function  DC 6/23, just a parfor package DC 5/23

p=gcp("nocreate");	%get the parpool number; 
lpar=p.NumWorkers;  	% how many seperate chromosomes to create 
[~,inC]=sort([C2.tfg]);	% should likely already be sorted 
C2=C2(inC);		% C2 input sorted
clear Cloc;		% not needed, but clarifies context for parfor 
 
parfor ipar=1:lpar          % use processors
    % use elite leaders
    il=1;

    for jj=ipar:lpar:numel(C2)
        Cloc(il)=C2(jj);	% Cloc's may not be all same size note 
        il=il+1;		% increment the Cloc 
    end

    [~, Cloc]=geneselect2(Ngen, Gmax,  Cloc, yERR, yOBS, CALC, maskCALC);	%the actual process 
    Cp(ipar,:)=Cloc; 			% has to be Cp(ipar,1:numel(Cloc))=Cloc; % non p
     
end

tmp=Cp(:);		% copy so that Cp untouched here 
tmp2={tmp.tfg};     	% list of tfg's (target functions)
mask=cellfun(@isempty,tmp2);  	%   any nulls?
tmp(mask)=[];			% remove them 

[~,inC]=sort([tmp.tfg]);    % sort the list, with no empties
Cfin=tmp(inC);              % return same size as input

return
end

function [tfg, C] = geneselect2 (Ngen,  Gmax, C, yERR, yOBS, CALC, maskCALC )
% execute Ngen generations for Gmax genes in chromosomes C with a target function 
% dependent on oberveds, yOBS, with variance yERR (reciprocal weights) compared to 
%  calculated values CALC -- a very large precalulated mx which may be partly masked by
%	maskCALC 
u_mask=false;   % using maskCALC or not? DC 6/23
comm=true; 	% put out comments in this case 
usegpu=isgpuarray(Ngen);	% test 
if nargin==7	% maskCALC provided 
    u_mask=true;
    if numel(maskCALC)~= size(CALC,1), error('geneselect2: bad mask size'); end
    XR=int32(find(~maskCALC));  % unsuppressed values
    % now randomize them 
    tmp=randperm(numel(XR)); 
    XR=int32(XR(tmp));
    comm=false;  % no comments lines with maskCALC 
end

function rchi = nchi (rr, yerr,y)
% normalze the raw chi to real chi  DC 7/23
yerr=yerr(:).^2;
y=y(:).^2;
rchi=rr.*sum(yerr)./sum(y)/(numel(yerr)-1); % provides a scale equivalent, parameterless number 
return
end

function [maskCALC , modC]= reduceCALC ( C, CALC,  maskCALC)  % modify to pool all
% filter out the below median set of contributors in CALC to chromosome tf into maskCALC 
% exclude values in C(1)  ; DC 4/23
[~,inC]=sort([C.tfg]);	
C=C(inC);   %probably presorted, but anyway...
tmpip=[C(:).ip];        % put the ip's in a regular array
Gmax=numel(C(1).ip);    % number of cols for ip
tmpwt=1./[C(:).tfg];    % weight of each contribution as tfg^-1
tmpwt=tmpwt./mean(tmpwt);   		% normalize to around 1.0
tmpip=reshape(tmpip, Gmax,numel(C) );   % MX form
Xpip=tmpip(:,1);        % first set
tmpM=size(CALC,1);
pool=zeros(tmpM,1,'single');    % vector for pool for each ip

for ii=1:numel(C)   % loop over C's
    w=tmpwt(ii);    % test weight
    ip_1=tmpip(:,ii);   	% test vector of ip's from C(i)
    maskw=pool(ip_1)< w;    	% any member < w?
    pool(ip_1(maskw))=w;    	% make up members to w value
end

% boost the Xpip set, so prime set not lost
pool(Xpip(:))=max(pool)+1; 	% lim could be 0
check=median(pool(pool>0));     % use the median, as mean could be skewed
mask = pool <= check & pool>0;  % what's less than median, but active
NewN=numel(find(~mask));
preN=numel(find(maskCALC));

if NewN  ==  preN
    modC=C;
    warning('reduceCALC 342 no reduction'); % should never be, but check
    return
end

frac=0.9;   % dont throw out everything
flim=check;
nzeros=nnz(pool);

while ( NewN)  < (1-frac)*(size(CALC,1)) % note only looking for exceptions in this group, modified for steps
    flim=(1-frac)*flim;
    mask=pool <= flim & pool > 0 ;
    NewN=numel(find(~mask)  );
    if NewN <=nzeros  % out of ranges
        flim=flim./(1+frac);
        mask = pool <=flim & pool>0;
        break;
    end
end
% increase if too few setout

target=min((size(CALC,1)-preN ).*frac, nzeros-Gmax-1);
range=[max(pool) , min(pool(pool>0))];
flim=check;
delf=1;

while nnz(mask) < target
    flim=range(2) + (flim-range(2))*(1+frac.*delf/10) ;
    mask=(pool <= flim) & pool > 0;
    delf=delf+frac;
end

delf=delf-frac;
flim=range(2) + (flim-range(2))*(1+frac.*delf/10) ;
mask=(pool <= flim) & pool > 0;

maskCALC=mask | maskCALC;

modC=C; % unchanged ... except for sort
return
end
 
tmp=[C(:).tfg];		% setp to sort the C set 
[~,inC]=sort(tmp);
C=C(inC);
 
% setup ranges
 
n_delta=floor(numel(C)./Gmax);  % set number of mutations involved 
oldtf=tmp;  			% carry a single vector
 
n_XR=(numel(XR));		% count of XR
n_CALC=(size(CALC,1));		% count of CALC
tfg=zeros(numel(C), Ngen,'single' );
[~,inC]=sort(oldtf);
OLDMIN=oldtf(inC(1));

for igen=1:Ngen  % outer loop for all generations 
    [~,inC]=sort(oldtf);	% sort withing the loop 
    C=C(inC);
    tmpip1=sort(int32(C(1).ip));  % for memory access improvement use the sort  

    if usegpu
        tmpip1=gpuArray(tmpip1);
    end

    ii=1;

    if oldtf(inC(1)) ~= OLDMIN   % improved this generation? 
        if comm, disp(num2str([igen inC(ii) oldtf(inC(1)) OLDMIN],'%5g ')); end
        OLDMIN=oldtf(inC(1));
    else
        n_delta=n_delta+1;   	% increase the loop size -- fewer max mutations 
    end

    for ii=2:numel(inC) % loop chromosomes 
        tmpip=tmpip1;   % copy the ip numbers 
        istage=min(max(1, floor(ii/n_delta)),Gmax);   	% how many mutations >= 1, <= Gmax 
        ipos=randperm(Gmax, istage );   		% position of mutations 

        if u_mask   % mask needs an indirect 
            tmp=randperm(n_XR,istage );	% pointer to locations to test 
             stage=XR(tmp);		% locations in CALC(maskCALC) 
        else
            stage=int32(randperm(n_CALC, istage));	% location in  CALC 
        end

        if numel(ipos) ~= numel(stage)			% debug test, can be removed 
            warning('bad stage'); disp(ipos);
            disp(istage);
            disp(stage);
        end

        tmpip(ipos)=stage;   		% substitute stage values in ipos 
        C(ii).ip=(int32(tmpip));  % update C stack 
    end

    tfg(:,igen)=gettf(yERR, C, yOBS, CALC); % calculate chi's 
    oldtf=tfg(:,igen);

    for jj=1:numel(C), C(jj).tfg=oldtf(jj); end   % update C values at end 

end
return
end

function FFd = FFdis(N , FFd )
%% test phe phe differences  DC 6 /23   can be used with 1 arg or cumulative
if nargin ==1 
    FFd=cell({});
else
    %FFd=FFd;  %nothing 
end
%ips=[1 100 1000 10000 100000 1000000]; close all;  % test line for code set
for kk=1:numel(N)	% over the N set 
    ii=N(kk);		% value sought
    if  ii <= numel(FFd) && numel(FFd{ii}), continue; end  % ensure in range and empty
    Ns=strip(int2str(ii)); % string version 

    if ispc  % compatible with PC or Linux environments 
        afil='Z:\Cowburn Lab\';
    else
        afil='~/shared/';
    end

  %  fil=['~/shared/TraDES-2/tmp12/FSFG12' Ns '_0000001.pdb'];
  	fil=[afil 'TraDES-2/tmp12/FSFG12' Ns '_0000001.pdb'];  % filename 
    d=dir(fil);	% check the directory 
    fil=[d.folder filesep d.name];	% formal full filename

    s=readlines(fil);		% string read 
    mask=startsWith(s,"ATOM");	% ATOM? 
    mask2=contains(s,'PHE');	% PHE ?
    mask3=endsWith(s,'H');	% H ? 
	if ~any(mask3), error('no hydrogens ', num2str(ii)); end 
    mask=mask & mask2 & mask3 ; % trifecta
    s=s(mask);
    writelines(s,'tmpphedi.pdb');	% temp file
    out=pdbread('tmpphedi.pdb');	% reread as structure 
  
    Atom=out.Model.Atom;
     
    co=(zeros(numel(Atom),3));
    for jj=1:numel(Atom)
        co(jj,1:3)=[Atom(jj).X Atom(jj).Y  Atom(jj).Z];	% get coordinates 
    end
  
    dr=hypot23u(co);	% get the upper mx of distances 
    dr=sqrt(dr);	% sqrt of same 
    FFd{ii}=single(dr); % put in cell.   Note cell used to minimize space as most values will never be used, and the addressing is efficient 
   
end

return
end

function Gout = phedis(N)%% test phe phe differences
% get the distribution of Phe-Phe atom distances from the precalculated pdbs 
%ips=[1 100 1000 10000 100000 1000000]; close all; % test case 
for ii=N
    Ns=strip(int2str(ii));	% to id pdb 
    fil=['../TraDES-2/tmp12/FSFG12' Ns '_0000001.pdb'];	% filename of pdb 
    s=readlines(fil);		% get the pdb 
    mask=startsWith(s,"ATOM");	% filter for ATOM lines 
    s=s(mask);
    writelines(s,'tmpphedi.pdb');	% write a tmp file for the pdbread function 
    %[out, file]=gettmppdb(ii);		% alternate quick/dirty approach 
    out=pdbread('tmpphedi.pdb');	% read the pdb file into the structure 
    delete('tmpphedi.pdb'); 		% cleanup 
   	Atom=out.Model.Atom;

  for jj=1:numel(Atom)	% could be completely vectorized but the file read/writes are slow 
      maskp(jj)=isequal(Atom(jj).resName,'PHE'); 	% PHE id 
      maksh(jj) = contains(Atom(jj).AtomName,'H');	% H atoms 
      co(jj,1:3)=[Atom(jj).X Atom(jj).Y  Atom(jj).Z];	% all coordinates
  end

 	mask=maskp & maksh;	% PHE & H
  	co=co(mask,:);	% select xyz
  	dr=zeros(size(co,1));	% distance mx
  	dr=hypot23u(co);	% call own routine or equivalent 
	f(ii)=figure;
	tmp=sqrt(dr(dr>0));	% sqrt of all values 
%tmp=tmp./sum(tmp);
	mtmp=max(tmp); 
	q(ii)=histogram(tmp,50, 'Normalization','probability');
	xlabel('D'); ylabel('Pop'); title ([ 'Phe H distances ' fil ]); 
	axis([0 150 0 0.07]);
end
Gout.f=f;	% filename 
Gout.q=q;  % for graphic debugging
return
end

function x = hypot23u(a, filt)
%--- x = hypot23(a);
% the distances between all cols of a are calculated
% note square is returned, since filetering may reduce need for sqrt
% calculation, upper portion only; designed to minimize rounding errors etc 
% filt may be used as a sparse filter as a maximum value
%--- dc 11/96, 11/08

i_a=size(a);    % sparse not useful generally 
if (i_a(2) ~= 3)
    error('hypot23: incorrect size of input');
end
x=zeros(i_a(1));
switch nargin   % no filter set 
    case 1
        for i=1:i_a(1) - 1
            z1=a(i,1);
            z2=a(i,2);
            z3=a(i,3);
            for j= i+1:i_a(1)
                x1=z1-a(j,1);
                xx=x1*x1;
                x1=z2-a(j,2);
                xx=xx+x1*x1;
                x1=z3-a(j,3);
                x(i,j)=xx+x1*x1;
            end
        end
    case 2
        for i =1:i_a(1) -1
            z1=a(i,1);
            z2=a(i,2);
            z3=a(i,3);
            for j = i+1:i_a(1)
                x1 = z1 - a(j,1);
                xx = x1*x1;
                if xx < filt
                    x1=z2-a(j,2);
                    xx=xx+x1*x1;
                    if xx < filt
                        x1=z3-a(j,3);
                        xx=xx+x1*x1;
                        if xx < filt
                            x(i,j) = xx;
                        end
                    end
                end
            end
        end
end


