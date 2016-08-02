function wmDrop_resampleReconstructions_rotateCoreg1(subj,VOIs,tpts_of_interest, hex_size)
%
%
% resamples all reconstrucitons within each condition with replacement
% n_iter times
%
root = load_root;
addpath(fullfile(root,'mFiles','gridfit'));

rand_seed = 2738382;
rng(rand_seed);

if nargin < 1
    subj = {'AP81','AP82','AP83','AI81','AI82','AI83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};
end

if nargin < 2 

    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS','SuperWMDrop'}; % FEFnew
end

if nargin < 3
    tpts_of_interest = {[3 4],[7 8]};
end

if nargin < 4
    hex_size = 7;
end

if ~iscell(tpts_of_interest)
    tpts_of_interest{1} = tpts_of_interest;
end

recon_str = 'trnAvg1'; 

n_resample_iter = 1000;


u_subj = unique(cellfun(@(s) s(1:end-1),subj,'uniformoutput',0));

use_first_half_180 = 1;





trials_per_superrun = 54;


nblank = length(subj) * length(VOIs) * 3 * trials_per_superrun * 12; % 3 superruns, 54 trials per superrun, 12 tpts

all_conds = nan(nblank,5);
all_tpts = nan(nblank,1);
all_subj = nan(nblank,1);
all_vois = nan(nblank,1);
all_180_half = nan(nblank,1);

% n_conds x n_VOIs
startidx = 1;
for ss = 1:length(subj)
    
    this_subj_id = find(strcmpi(u_subj,subj{ss}(1:end-1)));
    
    for vv = 1:length(VOIs)
        
        fn = sprintf('%swmDrop_recons/%s_%s_hex%i_%s_rotate_coreg1.mat',root,subj{ss},VOIs{vv},hex_size,recon_str);
        fprintf('loading %s...\n',fn);
        load(fn);
        
        if vv == 1 && ss == 1
            all_recons = nan(nblank,res_t*res_r);

            [gridt,gridr] = cart2pol(gridx,gridy);
        end
        
        thisidx = startidx:(startidx+size(recons_vec{1},1)-1);
        
        
        all_recons(thisidx,:) = recons_vec{1};

        all_conds(thisidx,:) = conds;
        all_tpts(thisidx) = tpts;
        all_subj(thisidx) = this_subj_id*ones(size(conds,1),1);
        all_vois(thisidx) = vv*ones(size(conds,1),1);
        
        % within each superrun (54 trials), identify first and second half
        % of 180 deg offset trials (all_conds(:,3)==3)
        
        
        n_super_runs = length(thisidx)/trials_per_superrun;
        whichhalf = nan(length(thisidx),1); % NaN for sep~=3, 1 for first half, 2 for second half
        
        for nn = 1:n_super_runs
            
            sridx = ((nn-1)*trials_per_superrun+1):(nn*trials_per_superrun);
            tmpconds = conds(sridx,:);
            for rr = 1:3
                tmp180 = find(tmpconds(:,3)==3 & tmpconds(:,2)==rr);
                whichhalf(sridx(tmp180(1:(length(tmp180)/2)))) = 1;
                whichhalf(sridx(tmp180( (1+(length(tmp180)/2)):length(tmp180)))) = 2;
            end

            clear sridx tmpconds tmp180;
        end
        
        all_180_half(thisidx) = whichhalf;
        
        
        
        startidx = thisidx(end)+1;
        clear whichhalf
    end
end

if use_first_half_180 == 1
    valididx = 1:(startidx-1);
    valididx = valididx(all_180_half(valididx)~=2);
else
    valididx = 1:(startidx-1);
end
    
all_recons = all_recons(valididx,:);
all_conds  = all_conds(valididx,:);
all_tpts = all_tpts(valididx);
all_subj = all_subj(valididx);
all_vois = all_vois(valididx);
all_180_half = all_180_half(valididx);

% average over polar angle arms (across ecc)
all_recons_tmp = reshape(all_recons.',[res_r,res_t,size(all_recons,1)]); % r x th x trials
all_recons_avgr = squeeze(mean(all_recons_tmp,1)).'; % trials x th
clear all_recons_tmp;

% label th
myt = reshape(gridt,res_r,res_t);
myt = myt(1,:);


% generate set of resampled trials (let's use the same set of indices for
% each condition, ROI, etc...)
n_trials = sum(all_conds(:,2)==1 & all_vois==1 & ismember(all_tpts,tpts_of_interest{1}));
resample_idx = nan(n_resample_iter,n_trials);
for ii = 1:n_resample_iter
    resample_idx(ii,:) = randsample(1:n_trials,n_trials,'true');
end

subj_str = sprintf('ALLn%i',length(u_subj));

%this_pool = parpool(8);

for tt = 1:length(tpts_of_interest)
    
    
    
    for vv = 1:length(VOIs)
        
        for cc = 1:3
            
            thisidx = find(all_conds(:,2)==cc & ismember(all_tpts,tpts_of_interest{tt}) & all_vois==vv);
            
            
            m_recon_vec = nan(n_resample_iter,size(all_recons_avgr,2));
            
            for ii = 1:n_resample_iter

                m_recon_vec(ii,:) = mean(all_recons_avgr(thisidx(resample_idx(ii,:)),:),1);
            end            
            
            tstr = [];
            for ii = 1:length(tpts_of_interest{tt})
                tstr = [tstr num2str(tpts_of_interest{tt}(ii))];
            end
            run_on = datestr(now,30);
            
            fn2s = sprintf('%swmDrop_recons/%s_%s_R%i_tpts%s_hex%i_Iter%i_%s_resampleRotateCoreg1.mat',root,subj_str,VOIs{vv},cc,tstr,hex_size,n_resample_iter,recon_str);
            fprintf('saving to %s...\n',fn2s);
            save(fn2s,'m_recon_vec','res_r','res_t','run_on','subj','resample_idx','rand_seed','myt');
            clear m_recon_vec thisidx tstr run_on ;
            
        
        end
    end
end




return

