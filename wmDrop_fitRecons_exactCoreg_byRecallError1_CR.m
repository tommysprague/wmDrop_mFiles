function wmDrop_fitRecons_exactCoreg_byRecallError1_CR(subj,VOIs,tpts_of_interest, hex_size)
%
% adapted from wmDrop_fitReconstructions_exactCoreg_byRecallError1_resample.m
%
% Loads data the same way as plotReconstructions of same suffix, but runs
% fits on each subject, each reconstruction using gridfit/finetune
%
% now we can/should fit all reconstructions simultaneously - so first need
% to compute a evalpts x n_recons matrix of all reconstructions we want to
% fit
%
% resamples all reconstrucitons within each condition with replacement
% n_iter times
%
% splits trials within each memory condition; subject based on recall error
% - simple median split - but one subj (AP, sess 2) has one more high-error
% than low-error trial due to odd number of tirals within each condition (3
% super-runs rather than 2) (also, all median splits done using trials used
% for analysis - so first half within each super-run for 180 deg separation
% trials)
%
% trnAvg is default now
%
% TS, 3/17/2015 - adapted from _resample.m

%addpath /usr/local/serenceslab/tommy/gridfit/;

if nargin < 1
    subj = {'AP81','AP82','AP83','AI81','AI82','AI83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};
    %subj = {'AP81','AP82','AP83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};
    %subj = {'AP82'};%,'AP82','AP83','AI81','AI82','AI83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83'};
end

if nargin < 2 
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS','SuperWMDrop'}; % FEFnew
    %VOIs = {'V1','V3A','IPS0','IPS2','sPCS','SuperWMDrop','V2','V3','V4','IPS1','IPS3'};
    %VOIs = {'sPCS','SuperWMDrop'};
    %VOIs = {'IPS1','IPS2','IPS3'}; % FEFnew
end

if nargin < 3
    %tpts_of_interest = {[3 4],[7 8]};
    tpts_of_interest = {[3 4],[7 8]};
end

if nargin < 4
    hex_size = 7;
end

if ~iscell(tpts_of_interest)
    tpts_of_interest{1} = tpts_of_interest;
end

root = load_root;%'/usr/local/serenceslab/tommy/wmDrop/';

n_resample_iter = 1000;

ecc_to_align = 3.5; % align to this point

myTR = 2.25;  % repetition time, in sec

mycolors = [0 127 66; 56 9 124; 28 69 68]./255;

u_subj = unique(cellfun(@(s) s(1:end-1),subj,'uniformoutput',0));

use_first_half_180 = 1;

condstr = {'R1','R2','R2d'};




trials_per_superrun = 54;

%tpts_of_interest = [2 3 4]; % to start with...

nblank = length(subj) * length(VOIs) * 3 * trials_per_superrun * 12; % 3 superruns, 54 trials per superrun, 12 tpts
% TODO: fix this...
if hex_size == 7
    %all_recons = nan(nblank,51^2);
    all_recons = nan(nblank,101^2);
else
    all_recons = nan(nblank,101^2);
end
all_conds = nan(nblank,5);
all_recallErr = nan(nblank,1);
all_recallMed = nan(nblank,1);
all_tpts = nan(nblank,1);
all_subj = nan(nblank,1);
all_vois = nan(nblank,1);
all_180_half = nan(nblank,1);

% n_conds x n_VOIs
startidx = 1;
for ss = 1:length(subj)
    
    this_subj_id = find(strcmpi(u_subj,subj{ss}(1:end-1)));
    
    fnbehav = sprintf('%swmDrop_combinedBehav/%s_wmDrop_combinedBehav.mat',root,subj{ss});
    fprintf('loading %s...\n',fnbehav);
    load(fnbehav);
    
    % divide into first/second half of super-run for 180 deg separation
    % trials
    
    n_super_runs = size(conditions,1)/trials_per_superrun;
    whichhalf = nan(size(conditions,1),1); % NaN for sep~=3, 1 for first half, 2 for second half
    
    for nn = 1:n_super_runs
        
        sridx = ((nn-1)*trials_per_superrun+1):(nn*trials_per_superrun);
        tmpconds = conditions(sridx,:);
        for rr = 1:3
            tmp180 = find(tmpconds(:,3)==3 & tmpconds(:,2)==rr);
            whichhalf(sridx(tmp180(1:(length(tmp180)/2)))) = 1;
            whichhalf(sridx(tmp180( (1+(length(tmp180)/2)):length(tmp180)))) = 2;
        end
        %all_sridx = thisidx(sridx); % this is where to put the values in
        clear sridx tmpconds tmp180;
    end
    
    
    % compute median split across all included trials (so only first half
    % of 180 deg separation distance trials)
    
    this_med_split = nan(size(conditions,1),1);
    
    % should compute median using only first half of 180 deg separation
    % trials (conds(:,3)==3)
    for cc = 1:3

        thisidx = conditions(:,2)==cc & whichhalf~=2;
        
        
        this_median = median(recall_error(thisidx));

        if mod(sum(thisidx),2)~=0
            fprintf('subj %s has odd number of trials, "bad error" trials oversampled slightly\n',subj{ss});
        end

        
        this_med_split(thisidx) = (recall_error(thisidx)>=this_median)+1;
        clear thisidx this_median;
        
    end
    
    
    clear conditions;
    
    for vv = 1:length(VOIs)
        
        fn = sprintf('%swmDrop_recons/%s_%s_hex%i_trnAvg_exact_coreg1.mat',root,subj{ss},VOIs{vv},hex_size);
        fprintf('loading %s...\n',fn);
        load(fn);
        
        
        thisidx = startidx:(startidx+size(recons_vec{1},1)-1);
        
        
        all_recons(thisidx,:) = recons_vec{1};
        all_conds(thisidx,:) = conds;
        all_recallErr(thisidx) = repmat(recall_error',length(unique(tpts)),1);
        
        all_recallMed(thisidx) = repmat(this_med_split,length(unique(tpts)),1);%(all_recallErr(thisidx)>this_median)+1;
        all_tpts(thisidx) = tpts;
        all_subj(thisidx) = this_subj_id*ones(size(conds,1),1);
        all_vois(thisidx) = vv*ones(size(conds,1),1);
        
        % within each superrun (54 trials), identify first and second half
        % of 180 deg offset trials (all_conds(:,3)==3)
        
        
%         n_super_runs = length(thisidx)/trials_per_superrun;
%         whichhalf = nan(length(thisidx),1); % NaN for sep~=3, 1 for first half, 2 for second half
        
%         for nn = 1:n_super_runs
%             
%             sridx = ((nn-1)*trials_per_superrun+1):(nn*trials_per_superrun);
%             tmpconds = conds(sridx,:);
%             for rr = 1:3
%                 tmp180 = find(tmpconds(:,3)==3 & tmpconds(:,2)==rr);
%                 whichhalf(sridx(tmp180(1:(length(tmp180)/2)))) = 1;
%                 whichhalf(sridx(tmp180( (1+(length(tmp180)/2)):length(tmp180)))) = 2;
%             end
%             %all_sridx = thisidx(sridx); % this is where to put the values in
%             clear sridx tmpconds tmp180;
%         end

        
        all_180_half(thisidx) = repmat(whichhalf,length(unique(tpts)),1);
        
        
        
        startidx = thisidx(end)+1;
        %clear whichhalf
    end
    clear whichhalf this_med_split;
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
all_recallErr = all_recallErr(valididx);
all_recallMed = all_recallMed(valididx);

%% ----- all fit parameters moved here ----

res = sqrt(size(all_recons,2));


maxecc = 6; % dva from fixation
%res = 101; % in x, y

% x, y, size (rad units), amp, base
ft_constr = [     -1 -3          0.25  -5 -5;
              maxecc  3 10/rad2fwhm(1)  10 10];


% --- fit parameters ---
%
% for now - using a square grid over full reconstruction, but can push this
% towards center later on
xydensity = 5; % 3 points per DVA
%xgrid = linspace(-maxecc,maxecc,xydensity*2*maxecc+1);
%ygrid = linspace(-maxecc,maxecc,xydensity*2*maxecc+1);

% match xgrid, ygrid to ft_constr above
%xgrid = linspace( ft_constr(1,1),ft_constr(2,1),xydensity*(ft_constr(2,1)-ft_constr(1,1))+1   );
%ygrid = linspace( ft_constr(1,2),ft_constr(2,2),xydensity*(ft_constr(2,2)-ft_constr(1,2))+1   );



% in rad units (so go up to FWHM = maxecc)
%sgrid = (0.5:0.25:maxecc)/rad2fwhm(1);
%sgrid = (0.5:0.25:maxecc)/rad2fwhm(1);
sgrid = ft_constr(1,3):0.25:ft_constr(2,3);


[gridx,gridy] = meshgrid(linspace(-maxecc,maxecc,res),linspace(-maxecc,maxecc,res));
gridx = reshape(gridx,numel(gridx),1);gridy = reshape(gridy,numel(gridy),1);
evalpts = [gridx gridy];

grid_params =  sgrid';

% need to move the below to the local maximum on every resampling iteration
% mygrid = make_grid(@make2dcos_grid,evalpts,grid_params);


pixwidth = maxecc*2/res;
ax = [];




% generate set of resampled trials (let's use the same set of indices for
% each condition, ROI, etc...)
n_trials = sum(all_conds(:,2)==1 & all_vois==1 & ismember(all_tpts,tpts_of_interest{1}) & all_recallMed == 1);
resample_idx = nan(n_resample_iter,n_trials);
for ii = 1:n_resample_iter
    resample_idx(ii,:) = randsample(1:n_trials,n_trials,'true');
end

subj_str = sprintf('ALLn%i',length(u_subj));

%this_pool = parpool(12);

for tt = 1:length(tpts_of_interest)
    
    
    
    for vv = 1:length(VOIs)
        
        for cc = 1:3
            
            for ms = [1 2] % loop over median splits (1 = low-error, 2 = high-error)
            
            thisidx = find(all_conds(:,2)==cc & ismember(all_tpts,tpts_of_interest{tt}) & all_vois==vv & all_recallMed==ms) ;
            
            ridx = 1;
            recons_to_fit = nan(size(evalpts,1),n_resample_iter);
            
            bf_grid = nan(n_resample_iter,size(ft_constr,2)+1);
            bf_fine = nan(n_resample_iter,size(ft_constr,2)+1);
            bffcn_grid = nan(n_resample_iter,size(evalpts,1));
            bffcn_fine = nan(n_resample_iter,size(evalpts,1));
            ex_flag_fine = nan(n_resample_iter,1);
            
            
            for ii = 1:n_resample_iter
                tic
                recons_to_fit(:,ridx) = mean(all_recons(thisidx(resample_idx(ii,:)),:),1)';
                
                
                % get local maximum
                [thismax,thismaxi] = max(recons_to_fit(:,ridx));
                thismaxX = gridx(thismaxi); thismaxY = gridy(thismaxi);
                %grid_params = repmat();
                
                mygrid = make_grid(@make2dcos_grid,evalpts,[thismaxX*ones(size(grid_params,1),1) thismaxY*ones(size(grid_params,1),1) grid_params]);
                
                fprintf('Resample iteration %i\n',ii);
                toc;
                [bf, err, bffcn] = gridfit(recons_to_fit(:,ii),mygrid,grid_params,[0 inf],0,0);
                tic;
                this_ft_constr = ft_constr;  % adjust constraints to keep x, y position within a single reconstruction pixel
                this_ft_constr(1,1) = thismaxX-pixwidth;
                this_ft_constr(1,2) = thismaxY-pixwidth;
                this_ft_constr(2,1) = thismaxX+pixwidth;
                this_ft_constr(2,2) = thismaxY+pixwidth;
                bf = [thismaxX*ones(size(bf,1),1) thismaxY*ones(size(bf,1),1) bf];
                
                [bf_ft, err_ft, bffcn_ft,ex_flag] = gridfit_finetune_constr(recons_to_fit(:,ii),@make2dcos_grid,bf,evalpts,this_ft_constr,0);
                
                
                
                
                
                bf_grid(ii,:) = [bf err]; bffcn_grid(ii,:) = bffcn; bf_fine(ii,:) = [bf_ft err_ft]; bffcn_fine(ii,:) = bffcn_ft'; ex_flag_fine(ii) = ex_flag;

                
                
                r_cond(ridx) = cc;
                r_subj(ridx) = ss;
                r_tpts(ridx) = tt;
                r_voi(ridx)  = vv;
                r_iter(ridx) = ii;
                r_med(ridx)  = ms;
                
                ridx = ridx+1;
                toc;
            end
            
        
            % do fitting here (across all resmapled iterations
            % simultaneously)
            %
            % then we don't need to re-sort the best fits later
            %[bf, err, bffcn] = gridfit(recons_to_fit,mygrid,grid_params,[0 inf]);
            %[bf_ft, err_ft, bffcn_ft,ex_flag] = gridfit_finetune_constr(recons_to_fit,@make2dcos_grid,bf,evalpts,ft_constr);

            
            % save best fits
            
            m_recon_vec = recons_to_fit';
            %bf_grid = [bf err']; bffcn_grid = bffcn; bf_fine = [bf_ft err_ft]; bffcn_fine = bffcn_ft'; ex_flag_fine = ex_flag;
            
            tstr = [];
            for ii = 1:length(tpts_of_interest{tt})
                tstr = [tstr num2str(tpts_of_interest{tt}(ii))];
            end
            run_on = datestr(now,30);
            
            fn2s = sprintf('%swmDrop_fits/%s_%s_R%i_tpts%s_ms%i_hex%i_Iter%i_fitExactCoreg_byRecallError1_CR.mat',root,subj_str,VOIs{vv},cc,tstr,ms,hex_size,n_resample_iter);
            fprintf('saving to %s...\n',fn2s);
            save(fn2s,'m_recon_vec','bf_grid','bffcn_grid','bf_fine','bffcn_fine','ex_flag_fine','sgrid','res','maxecc','run_on','ft_constr','subj');
            clear m_recon_vec thisidx bf_grid bffcn_grid bf_fine bffcn_fine ex_flag_fine tstr run_on;
            
            end
        
        end
    end
end



% now find each subj, ROI, etc and save out
%delete(this_pool);
% 
% for ss = 1:length(u_subj)
%     for vv = 1:length(VOIs)
%         for tt = 1:length(tpts_of_interest)
%             
%             % save out mean recon for each condition, fit after gridfit and
%             % fine-tune search, error, best-fit function for each condition
%             
%             thisidx = r_subj==ss & r_tpts == tt & r_voi == vv;
%             
%             m_recon_vec = recons_to_fit(:,thisidx)';
%             
%             cond = r_cond(thisidx);
%             
%             bf_grid = [bf(thisidx,:) err(thisidx)'];
%             bffcn_grid = bffcn(thisidx,:);
%             
%             bf_fine = [bf_ft(thisidx,:) err_ft(thisidx)];
%             bffcn_fine = bffcn_ft(:,thisidx)';
%             ex_flag_fine = ex_flag(thisidx);
%             tstr = [];
%             for ii = 1:length(tpts_of_interest{tt})
%                 tstr = [tstr num2str(tpts_of_interest{tt}(ii))];
%             end
%             
%             run_on = datestr(now,30);
%             
%             fn2s = sprintf('%swmDrop_fits/%s_%s_tpts%s_hex%i_fitExactCoreg1.mat',root,u_subj{ss},VOIs{vv},tstr,hex_size);
%             fprintf('saving to %s...\n',fn2s);
%             save(fn2s,'cond','m_recon_vec','bf_grid','bffcn_grid','bf_fine','bffcn_fine','ex_flag_fine','xgrid','ygrid','sgrid','res','maxecc','run_on','ft_constr');
%             clear m_recon_vec thisidx cond bf_grid bffcn_grid bf_fine bffcn_fine ex_flag_fine tstr run_on;
%         end
%     end
% end


return

