function wmDrop_plotTargetResp_resample1(subj,VOIs,tpts_of_interest,delay_tpts,hex_size)
%
% - extracts signals near the correct target position in PT and NPT
%   reconstructions (computed w/ computeReconstructions_exactCoreg1.m) and
%   plots them through time (like a HRF), as in Figure S4 of Sprague, Ester
%   & Serences, 2016
%
% NOTE: must include subj in order below to recover exact statistics
%
% First figure does not appear in SES, 2016 - this shows the time course of
% each target response individually. 
%
rand_seed = 21343;
rng(rand_seed);
save_stats = 0;
if nargin < 1
    subj = {'AP81','AP82','AP83','AI81','AI82','AI83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};
end

if nargin < 2 
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS','SuperWMDrop'};
end

if nargin < 3
    tpts_of_interest = 0:11;
end

if nargin < 4
    delay_tpts = {[3 4],[7 8]};
end

if nargin < 5
    hex_size = 7;
end

n_resample_iter = 1000;

root = load_root;

targ_pos = [3.5 0]; % align to this point
targ_window = 0.5; % extract data within this many dva of target position

myTR = 2.25;  % repetition time, in sec

load(fullfile(root,'wmDrop_colors.mat'));

u_subj = unique(cellfun(@(s) s(1:end-1),subj,'uniformoutput',0));

use_first_half_180 = 1;

load(fullfile(root,'wmDrop_subjSym.mat'));


condstr = {'R1','R2-neutral','R2-valid'};

maxecc = 6; % dva from fixation
res = 101; % in x, y

[gridx,gridy] = meshgrid(linspace(-maxecc,maxecc,res),linspace(-maxecc,maxecc,res));
gridx = reshape(gridx,numel(gridx),1);gridy = reshape(gridy,numel(gridy),1);

trials_per_superrun = 54;

%tpts_of_interest = [2 3 4]; % to start with...

nblank = length(subj) * length(VOIs) * 3 * 54 * length(tpts_of_interest); % 3 superruns, 54 trials per superrun, 12 tpts

all_conds = nan(nblank,5);
all_tpts = nan(nblank,1);
all_subj = nan(nblank,1);
all_vois = nan(nblank,1);
all_180_half = nan(nblank,1);
all_targ_resp = nan(nblank,2); % PT & NPT

% n_conds x n_VOIs
startidx = 1;
for ss = 1:length(subj)
    
    this_subj_id = find(strcmpi(u_subj,subj{ss}(1:end-1)));
    
    for vv = 1:length(VOIs)
        
        fn = sprintf('%swmDrop_recons/%s_%s_hex%i_trnAvg_exact_coreg1.mat',root,subj{ss},VOIs{vv},hex_size);
        load(fn);
        
        thisidx = startidx:(startidx+size(recons_vec{1},1)-1);
        
        if startidx == 1
            %res = sqrt(size(recons_vec{1},2));
            all_recons = nan(nblank,res^2);
            [gridx,gridy] = meshgrid(linspace(-maxecc,maxecc,res),linspace(-maxecc,maxecc,res));
            gridx = reshape(gridx,numel(gridx),1);gridy = reshape(gridy,numel(gridy),1);
            targ_idx = sqrt((gridx-targ_pos(1)).^2 + (gridy-targ_pos(2)).^2) < targ_window;
        end
        
        all_recons(thisidx,:) = recons_vec{1};
        all_conds(thisidx,:) = conds;
        all_tpts(thisidx) = tpts;
        all_subj(thisidx) = this_subj_id*ones(size(conds,1),1);
        all_vois(thisidx) = vv*ones(size(conds,1),1);
                
        all_targ_resp(thisidx,1) = mean(recons_vec{1}(:,targ_idx'),2);
        all_targ_resp(thisidx,2) = mean(recons_vec{2}(:,targ_idx'),2);
        % within each superrun (54 trials), identify first and second half
        % of 180 deg offset trials (all_conds(:,3)==3)
        
        
        n_super_runs = length(thisidx)/trials_per_superrun;
        whichhalf = nan(length(thisidx),1); % NaN for sep~=3, 1 for first half, 2 for second half
        
        
        for nn = 1:n_super_runs
            
            sridx = ((nn-1)*trials_per_superrun+1):(nn*trials_per_superrun);
            tmpconds = conds(sridx,:);
            for rr = 1:3
                tmp180 = find(tmpconds(:,3)==3&tmpconds(:,2)==rr);
                whichhalf(sridx(tmp180(1:(length(tmp180)/2)))) = 1;
                whichhalf(sridx(tmp180( (1+(length(tmp180)/2)):length(tmp180)))) = 2;
               
            end
            clear sridx tmpconds tmp180;
        end

        
        
        all_180_half(thisidx) = whichhalf;
        

        clear whichhalf
        
        startidx = thisidx(end)+1;
    end
end

valididx = 1:(startidx-1);
all_recons = all_recons(valididx,:);
all_conds  = all_conds(valididx,:);
all_tpts = all_tpts(valididx);
all_subj = all_subj(valididx);
all_vois = all_vois(valididx);
all_180_half = all_180_half(valididx);
all_targ_resp = all_targ_resp(valididx,:);



% 1 figure, 3 columns, each for R1/R2-neutral/R2-valid
% n_vois rows
% in each subplot,lines for each target (PT = solid, NPT = dashed)


ax = [];
figure;

for cc = 1:3   % column
    for vv = 1:length(VOIs)   % row
        
        
        ax(end+1) = subplot(length(VOIs),length(condstr),cc+(vv-1)*length(condstr));hold on;
        
        this_timeseries = nan(length(tpts_of_interest),2,n_resample_iter);
        

        n_trials = sum(all_conds(:,2)==cc & all_vois==vv & all_tpts==0 & (isnan(all_180_half)|(all_180_half==1)));
        

        this_data = nan(length(tpts_of_interest),2,n_trials);
        for tt = 1:length(tpts_of_interest)
            this_idx = all_tpts == tpts_of_interest(tt) & all_conds(:,2)==cc & all_vois==vv & (isnan(all_180_half)|(all_180_half==1));
            this_data(tt,:,:) = all_targ_resp(this_idx,:)';
        end
        
        % generate resmapling idx       
        resample_idx = nan(n_resample_iter,n_trials);
        for ii = 1:n_resample_iter
            resample_idx(ii,:) = randsample(1:n_trials,n_trials,'true');
        end
        
        
        
        for ii = 1:n_resample_iter           
            this_timeseries(:,:,ii) = mean(this_data(:,:,resample_idx(ii,:)),3);
        end
        
        mm = mean(this_timeseries,3);
        ee{1} = prctile(squeeze(this_timeseries(:,1,:)),[2.5 97.5]);
        ee{2} = prctile(squeeze(this_timeseries(:,2,:)),[2.5 97.5]);
        
        plot(tpts_of_interest*myTR, mm(:,1),'-' ,'color',mycolors(cc,:),'LineWidth',1.5);
        plot(tpts_of_interest*myTR, mm(:,2),'--','color',mycolors(cc,:),'LineWidth',1.5);
        
        
        if cc == 1
            ylabel(VOIs{vv});
        end
        
        if vv == 1
            title(condstr{cc});
        end
        
        if vv == length(VOIs) && cc == 1
            xlabel('Time (s)');
        end
        
    end
end


set(ax,'YLim',[-0.12 0.45],'Box','off','TickLength',[0.025 0.025],'XLim',myTR*[tpts_of_interest(1) tpts_of_interest(end)],'XTick',[0 8 16 24]);

%% plot a difference timeseries
% Fig. S4A-C of Sprague, Ester & Serences, 2016
% n_ROIs x n_conds, solid lines are mean over resampling iterations, dashed
% lines are 95% CI over resampling iterations
%
% NOTE: resampled a second time - to preserve stats as reported, must do
% the multiple resampling...sorry about that.
ax = [];
figure;

all_timeseries = cell(3,length(VOIs));
p_val_timeseries = cell(3,length(VOIs));

pidx = 1;
for cc = 1:3   % column
    for vv = 1:length(VOIs)   % row
                
        ax(end+1) = subplot(length(VOIs),length(condstr),cc+(vv-1)*length(condstr));hold on;
        
        this_timeseries = nan(length(tpts_of_interest),2,n_resample_iter);
        

        n_trials = sum(all_conds(:,2)==cc & all_vois==vv & all_tpts==0 & (isnan(all_180_half)|(all_180_half==1)));
        

        this_data = nan(length(tpts_of_interest),2,n_trials);
        for tt = 1:length(tpts_of_interest)
            this_idx = all_tpts == tpts_of_interest(tt) & all_conds(:,2)==cc & all_vois==vv & (isnan(all_180_half)|(all_180_half==1));
            this_data(tt,:,:) = all_targ_resp(this_idx,:)';
        end
        
        % generate resmapling idx
        resample_idx = nan(n_resample_iter,n_trials);
        for ii = 1:n_resample_iter
            resample_idx(ii,:) = randsample(1:n_trials,n_trials,'true');
        end

        
        for ii = 1:n_resample_iter
            
            this_timeseries(:,:,ii) = mean(this_data(:,:,resample_idx(ii,:)),3);
            
        end
        
        dd = squeeze(this_timeseries(:,1,:) - this_timeseries(:,2,:));
        mm = mean(dd,2);
        

        ee = prctile(dd',[2.5 97.5]);
        
        % leaving these as one-tailed, because we're testing for PT >> NPT
        p_val_timeseries{cc,vv} = mean(dd<=0,2);
        
        % also compute p-values for delay periods
        for tt = 1:length(delay_tpts)
            dd_delay_tmp(tt,:) = mean(dd(delay_tpts{tt},:),1);
            p_val_delay(pidx,tt) = mean(dd_delay_tmp(tt,:)<=0);
        end
        tmpdiff = dd_delay_tmp(2,:)-dd_delay_tmp(1,:);
        p_val_delay_comp(cc,vv) = 2*min(mean(tmpdiff>=0),mean(tmpdiff<=0));
        
        p_val_delay_labels(pidx,:) = [cc vv];
        
        pidx = pidx+1;
 
        clear dd_delay_tmp tmpdiff;
        

        
        plot(tpts_of_interest*myTR, mm,'-','color',mycolors(cc,:),'LineWidth',2.5);
        plot(tpts_of_interest*myTR, ee(1,:)','--','color',mycolors(cc,:),'LineWidth',1);
        plot(tpts_of_interest*myTR, ee(2,:)','--','color',mycolors(cc,:),'LineWidth',1);
        plot(myTR*[tpts_of_interest(1) tpts_of_interest(end)],[0 0],'k--');
        
        plot(8,-0.1,'k^','MarkerFaceColor','k','MarkerSize',5);
        
        
        
        if cc == 1
            ylabel(VOIs{vv});
        end
        
        if vv == 1
            title(condstr{cc});
        end
        
        if vv == length(VOIs) && cc == 1
            xlabel('Time (s)');
        end
        
        if vv~=length(VOIs)
            set(gca,'XTickLabel',[]);
            
        end
        if cc~=1
            set(gca,'YTickLabel',[]);
        end
        
        all_timeseries{cc,vv} = this_timeseries; % save this for stats below
        
    end
end

all_ylim = cell2mat(get(ax,'YLim'));

set(ax,'YLim',[-0.1 0.25],'Box','off','TickLength',[0.02 0.02],'XLim',myTR*[tpts_of_interest(1) tpts_of_interest(end)],'XTick',[0 8 16 24]);






%% also show a bar plot of each delay period


% 1 x n_ROIs
% each ROI has 3 bar gropus, 2 bars per group
figure;
axm = [];
for vv = 1:length(VOIs)
    for tt = 1:length(delay_tpts)
        axm(end+1) = subplot(length(VOIs),length(delay_tpts),(vv-1)*length(delay_tpts)+tt);hold on;
        
        this_targresp = nan(3,2,length(u_subj));
        
        for cc = 1:3
            for ss = 1:length(u_subj)
                
                
                if use_first_half_180 == 1
                    thisidx = all_subj==ss & all_conds(:,2)==cc & ismember(all_tpts,delay_tpts{tt}) & all_vois==vv & (isnan(all_180_half)|(all_180_half==1));
                else
                    thisidx = all_subj==ss & all_conds(:,2)==cc & ismember(all_tpts,delay_tpts{tt}) & all_vois==vv;
                end
                
                this_targresp(cc,1,ss) = mean(all_targ_resp(thisidx,1));
                this_targresp(cc,2,ss) = mean(all_targ_resp(thisidx,2));

                
                
                
            end
            
            bar([cc-.15],mean(this_targresp(cc,1,:)),'BarWidth',0.3,'EdgeColor',mycolors(cc,:),'FaceColor',mycolors(cc,:),'LineWidth',2);
            bar([cc+.15],mean(this_targresp(cc,2,:)),'BarWidth',0.3,'EdgeColor',mycolors(cc,:),'FaceColor',[1 1 1],'LineWidth',2);
            
            
            for ss = 1:length(u_subj)
                plot([cc-.15 cc+.15],squeeze(this_targresp(cc,:,ss)),'-','Color',[0.6 0.6 0.6],'LineWidth',1);
                plot(cc-.15,this_targresp(cc,1,ss),subjsym.(u_subj{ss}),'Color',mycolors(cc,:),'MarkerFaceColor',[1 1 1],'MarkerSize',6,'LineWidth',1);
                plot(cc+.15,this_targresp(cc,2,ss),subjsym.(u_subj{ss}),'Color',mycolors(cc,:),'MarkerFaceColor',[1 1 1],'MarkerSize',6,'LineWidth',1);
            end
            
            
        end
        
        if vv == 1
            title(sprintf('Delay %i',tt));
        end
        if tt == 1
            ylabel(VOIs{vv});
        end
        
        set(gca,'XTick',1:3,'XTickLabel',[]);
        if vv == length(VOIs)
            set(gca,'XTickLabel',condstr);
        end
        
        hold off;
    end
end

all_ylim = cell2mat(get(axm,'YLim'));

set(axm,'YLim',[-0.12 0.45],'Box','off','TickLength',[0.025 0.025],'XLim',[0.2 3.8]);


% FDR threshold - if 11 ROIs, we assume SuperWMDrop is one of them, so we
% exclude it from FDR computation (sorry, hacky...)
if length(VOIs)==11
    tmpp = p_val_delay(p_val_delay_labels(:,2)~=11,:);
    myfdr_delay(1) = fdr(tmpp(:),0.05); clear tmpp;
    tmpp = p_val_delay(p_val_delay_labels(:,2)==11,:);
    myfdr_delay(2) = fdr(tmpp(:),0.05); clear tmpp;
    
    tmpp = p_val_delay_comp(:,1:10);
    myfdr_delay_comp(1) = fdr(tmpp(:),0.05); clear tmpp;
    tmpp = p_val_delay_comp(:,11);
    myfdr_delay_comp(2) = fdr(tmpp(:),0.05); clear tmpp;

end

if save_stats == 1
    fn2s = sprintf('%swmDrop_stats/n%i_targResp_resample1_%s.mat',root,length(u_subj),datestr(now,30));
    fprintf('Saving stats to %s...\n',fn2s);
    
    save(fn2s,'p_val_timeseries','p_val_delay','p_val_delay_labels','p_val_delay_comp','myfdr_delay','myfdr_delay_comp','rand_seed','subj','VOIs','tpts_of_interest','delay_tpts','n_resample_iter');
end


return

