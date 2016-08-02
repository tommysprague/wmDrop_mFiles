function wmDrop_plotHRFs_ERA1(subj,VOIs,delay_timepoints)
%
% compute & plot HRFs (simple event-related average) and mean delay period activation
% (as plotted in Fig. S2 of Sprague, Ester & Serences, 2016)
% compare mean univariate response between delays (stats printed & saved)


if nargin < 1
    subj = {'AP81','AP82','AP83','AI81','AI82','AI83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};

end

if nargin < 2 
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS'};%,'iPCS','DLPFC','SMA'}; 

end
if nargin < 3
    delay_timepoints = {[3 4], [7 8]};
end


which_conds = [1 2 3];

root = load_root;

hrf_ylim = [-.3 0.5];
avg_ylim = [-0.3 0.3];

u_subj = unique(cellfun(@(s) s(1:end-1),subj,'uniformoutput',0));


myTR = 2.25;  % repetition time, in sec

load(fullfile(root,'wmDrop_colors.mat'));


%% extract all data

trials_per_superrun = 54;  % 54 when not noCue

nblank = trials_per_superrun*3*length(subj)*length(VOIs);


all_conds = nan(nblank,5);
all_subj = nan(nblank,1);
all_VOIs = nan(nblank,1);
all_180_half = nan(nblank,1);

startidx = 1;

for ss = 1:length(subj)
    this_subj_id = find(strcmpi(u_subj,subj{ss}(1:end-1)));
    for vv = 1:length(VOIs)
        fn = sprintf('%swmDrop_trialData/%s_%s_wmDrop_trialData.mat',root,subj{ss},VOIs{vv});
        fprintf('loading %s...\n',fn);
        thishrf = load(fn);
        
        if ss == 1 && vv == 1
            which_TRs = unique(thishrf.tpts);
            nTRs = length(which_TRs);
            
            all_mean_hrfs = nan(nblank,nTRs); % will mean over voxels to get these
        end
        
        % 1. mean across voxels
        % 2. reshape to n_trials x n_tpts
        %  all_mean_hrfs will be one hrf per trial (ERA)
        
        tmpm = mean(thishrf.tst,2);
        tmph = reshape(tmpm,size(thishrf.tst,1)/nTRs,nTRs);
        
        % because conds repeats over each time point, truncate to first
        tmpc = thishrf.conds(thishrf.tpts==which_TRs(1),:);
        
        thisidx = startidx:(startidx+size(tmph,1)-1);
        
        
        all_mean_hrfs(thisidx,:) = tmph;
        all_conds(thisidx,:) = tmpc;
        all_subj(thisidx) = this_subj_id;
        all_VOIs(thisidx) = vv;
        
        
        % within each superrun (54 trials), identify first and second half
        % of 180 deg offset trials (all_conds(:,3)==3)
        
        
        n_super_runs = length(thisidx)/trials_per_superrun;
        whichhalf = nan(length(thisidx),1); % NaN for sep~=3, 1 for first half, 2 for second half
        
        for nn = 1:n_super_runs
            
            sridx = ((nn-1)*trials_per_superrun+1):(nn*trials_per_superrun);
            tmpconds = tmpc(sridx,:);
            ru = unique(tmpconds);
            for rr = 1:length(ru) 
                tmp180 = find(tmpconds(:,3)==3&tmpconds(:,2)==ru(rr));
                whichhalf(sridx(tmp180(1:(length(tmp180)/2)))) = 1;
                whichhalf(sridx(tmp180( (1+(length(tmp180)/2)):length(tmp180)))) = 2;
                
            end
            clear sridx tmpconds tmp180;
        end
        
        all_180_half(thisidx) = whichhalf;
               
        
        startidx = thisidx(end)+1;
        clear thisidx thishrf whichhalf tmpm tmph tmpc;
    end
end



allidx = 1:(startidx-1);

all_mean_hrfs = all_mean_hrfs(allidx,:);
all_conds = all_conds(allidx,:);
all_subj = all_subj(allidx);
all_VOIs = all_VOIs(allidx);
all_180_half = all_180_half(allidx);



%% plot HRFs for each ROI
%
% want one row, one column for each ROI 




figure;
for vv = 1:length(VOIs)
    subplot(1,length(VOIs),vv);hold on;
    
    
    plot([0 (nTRs-1)*myTR],[0 0],'k-');
    
    for cc = 1:length(which_conds)
        this_hrf = nan(length(u_subj),nTRs);
        for ss = 1:length(u_subj)
            % mean over trials for each subject
            thisidx = all_conds(:,2)==which_conds(cc) & all_subj==ss & all_VOIs==vv;  % and all_180_half==1
            this_hrf(ss,:) = mean(all_mean_hrfs(thisidx,:),1);
        end
        
        % plot mean of HRF (for now just mean)
        plot([0:(nTRs-1)]*myTR,mean(this_hrf,1),'-','Color',mycolors(which_conds(cc),:));
        
    end
    
    if vv == 1
        ylabel('BOLD response (Z-score)');
        xlabel('Time (s)');
    end
    
    title(VOIs{vv});
    
    set(gca,'XTick',[0:4:nTRs*myTR],'TickLength',[0.025 0.025],'YTick',-1:0.2:1,'Box','off','XLim',[0 (nTRs-1)*myTR]);
    
    if vv~= 1
        set(gca,'YTickLabel',[]);
    end
    
    %hold off;
end

ax = get(gcf,'Children');
set(ax,'YLim',[hrf_ylim]);
%match_ylim(ax);

% add cue indicators
cue_time = 8.5; % hold/drop cue occurs here

for aa = 1:length(ax)
    cue_y = get(ax(aa),'YLim');
    cue_y = cue_y(1)+0.05;
    plot(ax(aa),cue_time,cue_y,'k^','MarkerFaceColor','k','MarkerSize',5);
end


%% plot mean of each delay ([3 4] or [7 8])
%
% before plotting, resample all trials with replacement and save confidence
% intervals (95%)


all_mean_delay = nan(size(all_mean_hrfs,1),length(delay_timepoints));

for tt = 1:length(delay_timepoints)
    all_mean_delay(:,tt) = mean(all_mean_hrfs(:,delay_timepoints{tt}),2);
end

n_resample_iter = 1000;
randseed = 237834738;
rng(randseed);

all_cis = nan(length(which_conds)*length(delay_timepoints)*length(VOIs),2);
all_ci_labels = nan(length(which_conds)*length(delay_timepoints)*length(VOIs),3); % label conditions, timepoint, ROI
all_resampled_data = nan(length(which_conds)*length(delay_timepoints)*length(VOIs),n_resample_iter);
myidx = 1;


for cc = 1:length(which_conds)

    % generate resample idx here, then use the same one for each set of
    % comparisons (so that same trials are always compared, jsut diffferent
    % data from same trials - from different timepoints or ROIs)
    n_trials_to_resample = sum(all_conds(:,2)==which_conds(cc) & all_VOIs==1);
    resample_idx{cc} = nan(n_resample_iter,n_trials_to_resample);
    for iter = 1:n_resample_iter
        resample_idx{cc}(iter,:) = randsample(n_trials_to_resample,n_trials_to_resample,true);
    end
    
    for tt = 1:length(delay_timepoints)
        for vv = 1:length(VOIs)
            
            thisidx = all_conds(:,2)==which_conds(cc) & all_VOIs==vv;
            thisdata = all_mean_delay(thisidx,tt);
            
            for iter = 1:n_resample_iter
                
                all_resampled_data(myidx,iter) = mean(thisdata(resample_idx{cc}(iter,:)));
                
            end
            
            all_cis(myidx,:) = prctile(all_resampled_data(myidx,:),[2.5 97.5]);
            all_ci_labels(myidx,:) = [cc,tt,vv];
            
            myidx = myidx+1;
        end
    end
end


xrange = 0.4;
dx = linspace(-xrange/2,xrange/2,length(which_conds));


% BELOW IS FOR MEAN & 95% CI 
figure;
for tt = 1:length(delay_timepoints)
    subplot(1,length(delay_timepoints),tt);
    hold on;
    plot([0 length(VOIs)+1],[0 0],'k-');
    for vv = 1:length(VOIs)
        for cc = 1:length(which_conds)
            
            mm = mean(all_resampled_data(all_ci_labels(:,1)==cc&all_ci_labels(:,2)==tt&all_ci_labels(:,3)==vv,:));
            thisci = all_cis(all_ci_labels(:,1)==cc&all_ci_labels(:,2)==tt&all_ci_labels(:,3)==vv,:);
            
            plot(vv+dx(cc)*ones(1,2),thisci,'-','color',mycolors(cc,:));

            plot(vv+dx(cc),mm,'o','color',mycolors(cc,:),'MarkerFaceColor',mycolors(cc,:));
            
        end
    end
    title(sprintf('Delay %i',tt));
    ylabel('BOLD response (Z-score)');
    set(gca,'XTick',1:length(VOIs),'XTickLabel',VOIs,'Box','off','XLim',[0.5 length(VOIs)+0.5]);
end

ax = get(gcf,'Children'); match_ylim(ax);

%% print out stats

cond_compare = [1 2; 2 3; 1 3];

allp = nan(length(delay_timepoints),length(VOIs),size(cond_compare,1));
allptab = nan(length(delay_timepoints)*length(VOIs)*size(cond_compare,1),4);
myidx = 1;
for tt = 1:length(delay_timepoints)
    for vv = 1:length(VOIs)
        for cc = 1:size(cond_compare,1)
            thisidx1 = all_ci_labels(:,1)==cond_compare(cc,1) & all_ci_labels(:,2)==tt & all_ci_labels(:,3)==vv;
            thisidx2 = all_ci_labels(:,1)==cond_compare(cc,2) & all_ci_labels(:,2)==tt & all_ci_labels(:,3)==vv;
            thisdiff = all_resampled_data(thisidx1,:) - all_resampled_data(thisidx2,:);
            thisp = 2*min(mean(thisdiff<0),mean(thisdiff>0));
            allp(tt,vv,cc)=thisp;
            allptab(myidx,:) = [tt vv cc thisp];
            myidx = myidx+1;
            fprintf('Delay %i\t%s\tR%i vs R%i:\tp = %0.03f\n',tt,VOIs{vv},cond_compare(cc,1),cond_compare(cc,2),thisp);
        end
        fprintf('\n');
    end
    fprintf('\n');
end

fdr_thresh = fdr(allp(:),0.05);
fprintf('Overall FDR threshold: %0.03f\n\n',fdr_thresh);

% save stats
fn2s = sprintf('%swmDrop_stats/n%i_meanDelayPeriod_resample_%iIter_%s.mat',root,length(u_subj),n_resample_iter,datestr(now,30));
fprintf('saving stats to %s...\n',fn2s);
save(fn2s,'allp','allptab','delay_timepoints','which_conds','VOIs','subj','u_subj','all_resampled_data','resample_idx','all_cis','all_ci_labels','fdr_thresh','randseed');


return