function wmDrop_fitReconstructions_exactCoreg1_CR(subj,VOIs,tpts_of_interest, hex_size)
%
% Compute resampled reconstructions and fits to resampled reconstructions
% as plotted in Fig. 7 of Sprague, Ester & Serences, 2016
%
% resamples all reconstrucitons within each condition with replacement
% n_iter times
%
%
% basically the same fitting procedure as CB 2014 - find local maximum on
% each fit iteration, search several possible sizes, each being GLM-fit. 
% 
%
% 
rng(787822979);
root = load_root;

% NOTE: to be numerically identical, must use subj names in order included
% below, and use 101x101 reconstructions (not computed by default). See
% wmDrop_allAnalyses1.m for further documentation. (and note that results
% will qualitatively align across parameters)
if nargin < 1
    subj = {'AP81','AP82','AP83','AI81','AI82','AI83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};
end

if nargin < 2 
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS','SuperWMDrop'}; 
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


n_resample_iter = 1000;



u_subj = unique(cellfun(@(s) s(1:end-1),subj,'uniformoutput',0));

use_first_half_180 = 1;



trials_per_superrun = 54;

%tpts_of_interest = [2 3 4]; % to start with...

nblank = length(subj) * length(VOIs) * 3 * trials_per_superrun * 12; % 3 superruns, 54 trials per superrun, 12 tpts

% if hex_size == 7
%     all_recons = nan(nblank,51^2);
% else
%     all_recons = nan(nblank,101^2);
% end
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
        
        fn = sprintf('%swmDrop_recons/%s_%s_hex%i_trnAvg_exact_coreg1.mat',root,subj{ss},VOIs{vv},hex_size);
        fprintf('loading %s...\n',fn);
        load(fn);
        
        
        if ss == 1 && vv == 1
            all_recons = nan(nblank,size(recons_vec{1},2));
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
            %all_sridx = thisidx(sridx); % this is where to put the values in
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

%% ----- all fit parameters moved here ----

res = sqrt(size(all_recons,2));


maxecc = 6; % dva from fixation

% GENERAL constraints for fine-tuning
% x, y, size (rad units), amp, base
ft_constr = [     NaN NaN          0.25  -5 -5;
                  NaN NaN 10/rad2fwhm(1)  10 10];
% above x,y are NaN because we impose iteration-specific constraints on x
% and y positions (allow a little wiggle; up to 1 pixel). those constraints
% are assigned below based on the best-fit gridfit.
              


% define size grid (span our fine-tuning constraints)
% in rad units (so go up to FWHM = maxecc)
sgrid = ft_constr(1,3):0.25:ft_constr(2,3);


[gridx,gridy] = meshgrid(linspace(-maxecc,maxecc,res),linspace(-maxecc,maxecc,res));
gridx = reshape(gridx,numel(gridx),1);gridy = reshape(gridy,numel(gridy),1);
evalpts = [gridx gridy];

grid_params =  sgrid';



pixwidth = maxecc*2/res;

% generate set of resampled trials (let's use the same set of indices for
% each condition, ROI, etc...)
n_trials = sum(all_conds(:,2)==1 & all_vois==1 & ismember(all_tpts,tpts_of_interest{1}));
resample_idx = nan(n_resample_iter,n_trials);
for ii = 1:n_resample_iter
    resample_idx(ii,:) = randsample(1:n_trials,n_trials,'true');
end

subj_str = sprintf('ALLn%i',length(u_subj));

%this_pool = parpool(8); with our heuristics, don't need to parpool -
%sufficient performance without

for tt = 1:length(tpts_of_interest)
    
        
    for vv = 1:length(VOIs)
        
        for cc = 1:3  % loop over conditions
            
            thisidx = find(all_conds(:,2)==cc & ismember(all_tpts,tpts_of_interest{tt}) & all_vois==vv);
            
            ridx = 1;
            recons_to_fit = nan(size(evalpts,1),n_resample_iter);
            
            bf_grid = nan(n_resample_iter,size(ft_constr,2)+1);
            bf_fine = nan(n_resample_iter,size(ft_constr,2)+1);
            bffcn_grid = nan(n_resample_iter,size(evalpts,1));
            bffcn_fine = nan(n_resample_iter,size(evalpts,1));
            ex_flag_fine = nan(n_resample_iter,1);
            
            for ii = 1:n_resample_iter
                recons_to_fit(:,ridx) = mean(all_recons(thisidx(resample_idx(ii,:)),:),1)';
                
                % get local maximum - this is used as 'fixed' center
                % coordinate for this grid search  through size. on each
                % resampling iteration, we linearly fit a*f(r)+b, where
                % f(r) is defined by each of a densely-sampled set of sizes
                % and maximum-defined position
                [thismax,thismaxi] = max(recons_to_fit(:,ridx));
                thismaxX = gridx(thismaxi); thismaxY = gridy(thismaxi);
                
                
                mygrid = make_grid(@make2dcos_grid,evalpts,[thismaxX*ones(size(grid_params,1),1) thismaxY*ones(size(grid_params,1),1) grid_params]);
                
                fprintf('Resample iteration %i\n',ii);
                [bf, err, bffcn] = gridfit(recons_to_fit(:,ii),mygrid,grid_params,[0 inf],0,0);
                
                % adjust constraints to keep x, y position within a single reconstruction pixel
                this_ft_constr = ft_constr;  
                this_ft_constr(1,1) = thismaxX-pixwidth;
                this_ft_constr(1,2) = thismaxY-pixwidth;
                this_ft_constr(2,1) = thismaxX+pixwidth;
                this_ft_constr(2,2) = thismaxY+pixwidth;
                bf = [thismaxX*ones(size(bf,1),1) thismaxY*ones(size(bf,1),1) bf];
                [bf_ft, err_ft, bffcn_ft,ex_flag] = gridfit_finetune_constr(recons_to_fit(:,ii),@make2dcos_grid,bf,evalpts,this_ft_constr,0);
            
                
                
                
                
                bf_grid(ii,:) = [bf err]; bffcn_grid(ii,:) = bffcn; bf_fine(ii,:) = [bf_ft err_ft]; bffcn_fine(ii,:) = bffcn_ft'; ex_flag_fine(ii) = ex_flag;
                
               
                ridx = ridx+1;
            end
            
            m_recon_vec = recons_to_fit';
            
            
            tstr = [];
            for ii = 1:length(tpts_of_interest{tt})
                tstr = [tstr num2str(tpts_of_interest{tt}(ii))];
            end
            run_on = datestr(now,30);
            
            fn2s = sprintf('%swmDrop_fits/%s_%s_R%i_tpts%s_hex%i_Iter%i_fitExactCoreg1_CR.mat',root,subj_str,VOIs{vv},cc,tstr,hex_size,n_resample_iter);
            fprintf('saving to %s...\n',fn2s);
            save(fn2s,'m_recon_vec','bf_grid','bffcn_grid','bf_fine','bffcn_fine','ex_flag_fine','sgrid','res','maxecc','run_on','ft_constr','subj');
            clear m_recon_vec thisidx bf_grid bffcn_grid bf_fine bffcn_fine ex_flag_fine tstr run_on;
            
        
        end
    end
end

try
    delete(this_pool);
end


return

