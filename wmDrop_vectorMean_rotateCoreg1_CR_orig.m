function wmDrop_vectorMean_rotateCoreg1_CR_orig(subj,VOIs,tpts_of_interest, hex_size)
%
% TCS 6/23/2015
% updated 8/15/2015 to plot comparison between represntational fidelity
% from delay 1 to delay 2, and save out stats
%
%for full time-course - can then group tpts, etc.
%
% computes a 1-d reconstruction along 3.5 dva ring (from rotateCoreg1.mat
% files)
%
%
% resamples all reconstrucitons within each condition with replacement
% n_iter times
%
% ts apr 14, 2016  - adding difference scores: (v1: delay 2-1) - v2: delay
% 2-1) - differences in "restoration" or "cue-related" effect
%
% ts apr 18, 2016 - added difference comparisons for each timepoint pair,
% stats (shows restoration effect is robust across delay periods chosen)
%

root = load_root;
addpath(fullfile(root,'mFiles','gridfit'));

% RNG seed used for Sprague, Ester & Serences, 2016: 2738382
rand_seed = 2738382;
rng(rand_seed);

% if "1" will save a stats file into wmDrop_stats
save_stats = 1;

if nargin < 1

    subj = {'AI81','AI82','AI83','AL81','AL82','AL83','AP81','AP82','AP83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};
    
end

if nargin < 2
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS','SuperWMDrop'}; % FEFnew
end

if nargin < 3
    tpts_of_interest = 0:9;
end

if nargin < 4
    hex_size = 7;
end

% let's also average vec_means at these timepoints and do stats for the 1-d
% reconstruction part of figure(s)
delay_tpts = {[3 4],[7 8]};




recon_str = 'trnAvg1'; 

n_resample_iter = 1000;



u_subj = unique(cellfun(@(s) s(1:end-1),subj,'uniformoutput',0));

use_first_half_180 = 1;

condstr = {'R1','R2-neutral','R2-valid'};




trials_per_superrun = 54;


nblank = length(subj) * length(VOIs) * 3 * trials_per_superrun * length(tpts_of_interest); % 3 superruns, 54 trials per superrun, 12 tpts

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
            %all_recons = nan(nblank,res_t);
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


all_recons_tmp = reshape(all_recons.',[res_r,res_t,size(all_recons,1)]); % r x th x trials
all_recons_avgr = squeeze(mean(all_recons_tmp,1)).'; % trials x th
clear all_recons_tmp;



myt = reshape(gridt,res_r,res_t);
myt = myt(1,:);
myx = (myt.'/(2*pi)) * 2*pi*mean(gridr); 


% generate set of resampled trials (let's use the same set of indices for
% each condition, ROI, etc...)
n_trials = sum(all_conds(:,2)==1 & all_vois==1 & all_tpts==tpts_of_interest(1));
resample_idx = nan(n_resample_iter,n_trials);
for ii = 1:n_resample_iter
    resample_idx(ii,:) = randsample(1:n_trials,n_trials,'true');
end

subj_str = sprintf('ALLn%i',length(u_subj));

all_vecmean = nan(length(VOIs)*3*length(tpts_of_interest),2,n_resample_iter);
all_labels = nan(length(VOIs)*3*length(tpts_of_interest),3);
%this_pool = parpool(8);
ridx = 1;
didx = 1;


all_vecmean_delay = nan(length(VOIs)*3*length(delay_tpts),2,n_resample_iter);
all_labels_delay = nan(length(VOIs)*3*length(delay_tpts),3);



for vv = 1:length(VOIs)
    
    for cc = 1:3
        for tt = 1:length(tpts_of_interest)
            
            thisidx = find(all_conds(:,2)==cc & all_tpts==tpts_of_interest(tt) & all_vois==vv);
            
            
            recons_to_fit = nan(size(all_recons_avgr,2),n_resample_iter);
            
            
            for ii = 1:n_resample_iter
                recons_to_fit(:,ii) = mean(all_recons_avgr(thisidx(resample_idx(ii,:)),:),1)';
            end
            
            
            
            m_recon_vec = recons_to_fit.';
            
            
            recon_vec_mean = [ mean( m_recon_vec.*cos(repmat(myt,n_resample_iter,1)),2) mean( m_recon_vec.*sin(repmat(myt,n_resample_iter,1)),2) ];
            
            all_vecmean(ridx,:,:) = recon_vec_mean.';
            all_labels(ridx,:) = [vv,cc,tt];
            
           
            ridx = ridx+1;
            
            
            
            clear m_recon_vec thisidx tstr run_on recon_vec_mean;
            
            
        end
        
        % do the same thing, but average over timepoints subtending each
        % delay
        for tt = 1:length(delay_tpts)
            thisidx = find(all_conds(:,2)==cc & ismember(all_tpts,delay_tpts{tt}) & all_vois==vv);
            recons_to_fit = nan(size(all_recons_avgr,2),n_resample_iter);
            
            
            for ii = 1:n_resample_iter
                recons_to_fit(:,ii) = mean(all_recons_avgr(thisidx(resample_idx(ii,:)),:),1)';
            end
            
            
            
            m_recon_vec = recons_to_fit.';
            
            
  
            recon_vec_mean = [ mean( m_recon_vec.*cos(repmat(myt,n_resample_iter,1)),2) mean( m_recon_vec.*sin(repmat(myt,n_resample_iter,1)),2) ];
            
            all_vecmean_delay(didx,:,:) = recon_vec_mean.';
            all_labels_delay(didx,:) = [vv,cc,tt];
            
           
            didx = didx+1;
            
            
            
            clear m_recon_vec thisidx tstr run_on recon_vec_mean;
            
        end
    end
end



%% do stats
allp = nan(length(VOIs),3,length(tpts_of_interest));
for vv = 1:length(VOIs)
    for cc = 1:3
        for tt = 1:length(tpts_of_interest)
            thisidx = all_labels(:,1)==vv & all_labels(:,2)==cc & all_labels(:,3)==tt;
            
            % ONE-TAILED!!!! doesn't make sense to look for -x
            allp(vv,cc,tt) = mean(squeeze(all_vecmean(thisidx,1,:))<0);
        end
    end
end

allp_delay = nan(length(VOIs),3,length(delay_tpts));
for vv = 1:length(VOIs)
    for cc = 1:3
        for tt = 1:length(delay_tpts)
            thisidx = all_labels_delay(:,1)==vv & all_labels_delay(:,2)==cc & all_labels_delay(:,3)==tt;
            
            % ONE-TAILED!!!! doesn't make sense to look for -x
            allp_delay(vv,cc,tt) = mean(squeeze(all_vecmean_delay(thisidx,1,:))<0);
        end
    end
end


% do stats comparing 2nd delay & first delay

allp_delay_compare = nan(length(VOIs),3); % compare delay 2 vs delay 1
for vv = 1:length(VOIs)
    for cc = 1:3
        %        for tt = 1:length(delay_tpts)
        thisidx1 = all_labels_delay(:,1)==vv & all_labels_delay(:,2)==cc & all_labels_delay(:,3)==1;
        thisidx2 = all_labels_delay(:,1)==vv & all_labels_delay(:,2)==cc & all_labels_delay(:,3)==2;
        
        % TWO-TAILED
        allp_delay_compare(vv,cc) = 2*min( mean(squeeze( all_vecmean_delay(thisidx1,1,:) < all_vecmean_delay(thisidx2,1,:)  )) , mean(squeeze( all_vecmean_delay(thisidx2,1,:) < all_vecmean_delay(thisidx1,1,:)  ))   );
        %        end
    end
end


load(fullfile(root,'wmDrop_colors.mat'));

my_ylim = [inf -inf];

myalpha = 0.05; % for trends

if length(VOIs)==11 % if including SuperWMDrop
    tmpp1 = allp(1:10,:,:);
    tmpp2 = allp(11,:,:);
    mythresh = [fdr(tmpp1(:),0.05)*ones(length(VOIs)-1,1); fdr(tmpp2(:),0.05)];  % for significance
    tmpp1_delay = allp_delay(1:10,:,:);
    tmpp2_delay = allp_delay(11,:,:);
    mythresh_delay = [fdr(tmpp1_delay(:),0.05)*ones(length(VOIs)-1,1); fdr(tmpp2_delay(:),0.05)];  % for significance
    
    tmpp1_delay_compare = allp_delay_compare(1:10,:);
    tmpp2_delay_compare = allp_delay_compare(11,:);
    mythresh_delay_compare = [fdr(tmpp1_delay_compare(:),0.05)*ones(length(VOIs)-1,1); fdr(tmpp2_delay_compare(:),0.05)];
    
    
else
    mythresh = fdr(allp(:),0.05)*ones(length(VOIs),1);
    mythresh_delay = fdr(allp_delay(:),0.05)*ones(length(VOIs),1);
    mythresh_delay_compare = fdr(allp_delay_compare(:),0.05)*ones(length(VOIs),1);
end

if save_stats == 1
    fn_stats = sprintf('%swmDrop_stats/n%i_vectorMean_thruTimeOrig_%sthru%s_%s_%s.mat',root,length(u_subj),VOIs{1},VOIs{end},recon_str,datestr(now,30));
    fprintf('saving stats to %s...\n',fn_stats);
    save(fn_stats,'mythresh','mythresh_delay','mythresh_delay_compare','allp','allp_delay','allp_delay_compare','tpts_of_interest','VOIs','condstr','n_resample_iter','rand_seed','delay_tpts');
end


%% plot vector mean (representational fidelity) timecourse
%
% Figure 5

myTR = 2.25;
figure;
for vv = 1:length(VOIs)
    thisax = [];
    my_ylim = [inf -inf];
    for cc = 1:3
        thisax(end+1) = subplot(length(VOIs),3,cc+3*(vv-1));hold on;
        
        % plot x axis
        plot(tpts_of_interest([1 end])*myTR,[0 0],'k-');
        
        this_tc = nan(length(tpts_of_interest),1);
        this_ci = nan(length(tpts_of_interest),2);
        for tt = 1:length(tpts_of_interest)
            thisidx = all_labels(:,1)==vv & all_labels(:,2)==cc & all_labels(:,3)==tt;
            this_tc(tt) = mean(squeeze(all_vecmean(thisidx,1,:)));
            this_ci(tt,:) = prctile(squeeze(all_vecmean(thisidx,1,:)),[2.5 97.5]);
            if allp(vv,cc,tt) <= mythresh(vv)   % filled markers
                plot(tpts_of_interest(tt)*myTR,0,'o','color',mycolors(cc,:),'MarkerFaceColor',mycolors(cc,:),'MarkerSize',3,'LineWidth',1);
            elseif allp(vv,cc,tt) <= myalpha % open markers
                plot(tpts_of_interest(tt)*myTR,0,'o','color',mycolors(cc,:),'MarkerSize',3,'MarkerFaceColor',[1 1 1],'LineWidth',1);
            end
            
        end
        
        plot(tpts_of_interest*myTR,this_tc,'-','Color',mycolors(cc,:),'LineWidth',1);
        plot(tpts_of_interest*myTR,this_ci.','--','Color',mycolors(cc,:),'LineWidth',0.5);
        
        if max(this_ci(:)) > my_ylim(2)
            my_ylim(2) = max(this_ci(:))+0.01;
        end
        
        if min(this_ci(:)) < my_ylim(1)
            my_ylim(1) = min(this_ci(:))-0.01;
        end
        
        if cc == 1
            ylabel(VOIs{vv});
        else
            set(gca,'YTickLabel',[]);
        end
        
        if vv == 1
            title(condstr{cc});
        end
        
        if vv~=length(VOIs)
            set(gca,'XTickLabel',[])
        end
        
        hold off;
    end
    %match_ylim(thisax);
    
    
    set(thisax,'YLim',my_ylim);
    set(thisax,'XLim',myTR*[tpts_of_interest(1) tpts_of_interest(end)],'XTick',tpts_of_interest(1):4:(myTR*tpts_of_interest(end)),'Box','off','TickLength',[0.025 0.025]);
end

set(gcf,'Position',[440    30   317   768]);

%% plot the vector mean (representational fidelity) for each delay
%
% Figure 6


figure;thisax = [];my_ylim = [inf -inf];xoffset = linspace(-0.15,0.15,2); % left/right this far
for cc = 1:3
    
    
    thisax(end+1) = subplot(3,1,cc);hold on;
    
    % plot an x axis
    plot([0 length(VOIs)+1],[0 0],'k-');
    
    for vv = 1:length(VOIs)
        this_vecmean = nan(length(delay_tpts),1);
        this_ci = nan(length(delay_tpts),2);
        
        for tt = 1:length(delay_tpts)
            thisidx = all_labels_delay(:,1)==vv & all_labels_delay(:,2)==cc & all_labels_delay(:,3)==tt;
            this_vecmean(tt) = mean(squeeze(all_vecmean_delay(thisidx,1,:)));
            this_ci(tt,:) = prctile(squeeze(all_vecmean_delay(thisidx,1,:)),[2.5 97.5]);
            
            % CI
            plot((vv+xoffset(tt))*[1 1],this_ci(tt,:),'-','LineWidth',1,'Color',mycolors(cc,:));
            
            if tt == 1
                plot(vv+xoffset(tt),this_vecmean(tt),'o','MarkerSize',7,'LineWidth',1,'Color',mycolors(cc,:),'MarkerFaceColor',mycolors(cc,:));
            else
                plot(vv+xoffset(tt),this_vecmean(tt),'o','MarkerSize',7,'LineWidth',1,'Color',mycolors(cc,:),'MarkerFaceColor',[1 1 1]);
            end
            
        end
        
    end
    
    set(gca,'XTick',1:length(VOIs),'XTickLabel',VOIs,'XLim',[0 length(VOIs)+1]);
    title(condstr{cc});
    if cc == 2
        ylabel('Representational fidelity (BOLD Z-score)');
    end
    
end

%match_ylim(thisax);
set(thisax,'YTick',[0 .05],'YLim',[-0.02 0.075]);


load(fullfile(root,'rwb.mat'));

%% compare between regions
% this works out to be a test for an interaction between delay and ROI for
% a pair of ROIs
%
% Figure S5

% can't do this analysis w/ 1 ROI...
if length(VOIs) > 1
    
    
% for this: exclude superWMDrop
    if strcmpi(VOIs{end},'SuperWMDrop') && length(VOIs)~=1
        oVOIs = VOIs;
        VOIs = {VOIs{1:(end-1)}};
    end
    
    
    
    
    thisp_diffm =cell(3,1);
    fdr_diffs = nan(3,1);
    figure;
    
    for cc = 1:3
        
        % I need a variable keeping track of differences between A:D2-D1 and B:
        % D2-D1
        
        this_diffm = nan(length(VOIs),length(VOIs));
        thisp_diffm{cc} = nan(length(VOIs),length(VOIs));
        
        for vv1 = 2:length(VOIs)
            
            % VOI 1, delay 1
            thisidxv1d1 = all_labels_delay(:,1)==vv1 & all_labels_delay(:,2)==cc & all_labels_delay(:,3)==1;
            % VOI 1, delay 2
            thisidxv1d2 = all_labels_delay(:,1)==vv1 & all_labels_delay(:,2)==cc & all_labels_delay(:,3)==2;
            
            % dist of their differences
            thisdiff_v1 = squeeze(all_vecmean_delay(thisidxv1d2,1,:) - all_vecmean_delay(thisidxv1d1,1,:) );
            % n_samp_iters long
            
            for vv2 = 1:(vv1-1)
                
                % VOI 1, delay 1
                thisidxv2d1 = all_labels_delay(:,1)==vv2 & all_labels_delay(:,2)==cc & all_labels_delay(:,3)==1;
                % VOI 1, delay 2
                thisidxv2d2 = all_labels_delay(:,1)==vv2 & all_labels_delay(:,2)==cc & all_labels_delay(:,3)==2;
                
                % dist of their differences
                thisdiff_v2 = squeeze(all_vecmean_delay(thisidxv2d2,1,:) - all_vecmean_delay(thisidxv2d1,1,:) );
                
                
                
                % do delay 2 - delay 1
                this_diffm(vv1,vv2) = mean(thisdiff_v2-thisdiff_v1);
                thisp_diffm{cc}(vv1,vv2) = 2*min(mean((thisdiff_v2-thisdiff_v1)<0),mean( (thisdiff_v2-thisdiff_v1)>0 ));
                
            end
        end
        
        diffax(cc) = subplot(1,3,cc);hold on;
        diff2plot = this_diffm; diff2plot(isnan(diff2plot))=0;
        imagesc(diff2plot);set(gca,'XTick',1:length(VOIs),'YTick',1:length(VOIs)','XTickLabel',VOIs,'YTickLabel',VOIs);
        colormap(rwb);
        xlabel('ROI 2');
        ylabel('ROI 1');
        title(sprintf('%s: ROI 2 - ROI 1',condstr{cc}));
        set(gca,'XTickLabelRotation',90)
        
        fdr_diffs(cc) = fdr(thisp_diffm{cc}(~isnan(thisp_diffm{cc})),0.05);
        
        fprintf('%s = FDR thresh = %0.03f\n',condstr{cc},fdr_diffs(cc));
        
        
        
        % draw square around relevant cells, black if FDR significant,
        % otherwise, (trend), gray square,
        for vv1 = 2:length(VOIs)
            
            for vv2 = 1:(vv1-1)
                
                if thisp_diffm{cc}(vv1,vv2)<=fdr_diffs(cc)
                    rectangle('Position',[(vv2-0.49) (vv1-0.49) 0.98 0.98],'EdgeColor',[0 0 0],'LineWidth',1);
                elseif thisp_diffm{cc}(vv1,vv2)<=0.05
                    rectangle('Position',[(vv2-0.49) (vv1-0.49) 0.98 0.98],'EdgeColor',[0.5 0.5 00.5],'LineWidth',1);
                end
                
            end
        end
        
        axis square ij;
        
    end
    
    
    myclim = cell2mat(get(diffax,'CLim'));
    newclim = max(abs(myclim(:)));
    set(diffax,'CLim',[-1 1]*newclim,'XLim',[0.5 length(VOIs)+0.5],'YLim',[0.5 length(VOIs)+0.5]);

    
    if save_stats == 1
        fn_stats = sprintf('%swmDrop_stats/n%i_vectorMean_diffsOrig_%sthru%s_%s_%s.mat',root,length(u_subj),VOIs{1},VOIs{end},recon_str,datestr(now,30));
        fprintf('saving stats to %s...\n',fn_stats);
        save(fn_stats,'fdr_diffs','thisp_diffm','VOIs','n_resample_iter','rand_seed','delay_tpts');%,'mythresh_delay_compare','allp','allp_delay','allp_delay_compare','tpts_of_interest','VOIs','condstr','n_resample_iter','rand_seed','delay_tpts');
    end
end

%% timepoint vs. timepoint comparisons, within each condition/ROI
%
% Figure S6
% 
% plotted as n_ROIs x n_conditions images, like above. FDR correction
% within-condition? perhaps withnin condition/ROI pair - huge number of
% comparisons here
all_difft = cell(length(VOIs),3); % these are the images we'll plot
allp_difft =cell(length(VOIs),3); % and corresponding p-values
fdr_difft = nan(3,1);
figure; axt = nan(length(VOIs),3);
for cc = 1:3
    for vv = 1:length(VOIs)
        
        
        
        all_difft{vv,cc} = nan(length(tpts_of_interest),length(tpts_of_interest));
        allp_difft{vv,cc} = nan(length(VOIs),length(VOIs));
        
        for tt2 = 2:length(tpts_of_interest)
            
            
            for tt1 = 1:(tt2-1)
                
                % first time point
                thisidxt1 = all_labels(:,1)==vv & all_labels(:,2)==cc & all_labels(:,3)==tt1;
                
                % paired with this timepoint
                thisidxt2 = all_labels(:,1)==vv & all_labels(:,2)==cc & all_labels(:,3)==tt2;
                
                
                % dist of their differences (tpt 2 - tpt 1)
                thisdifft = squeeze(all_vecmean(thisidxt2,1,:) - all_vecmean(thisidxt1,1,:) );
                % n_samp_iters long
                
                % do tpt 2 - tpt 1
                all_difft{vv,cc}(tt2,tt1) = mean(thisdifft);
                allp_difft{vv,cc}(tt2,tt1) = 2*min(mean(thisdifft<0),mean(thisdifft>0));
                
                clear thisidxt2 thisidxt1 thisdifft;
                
            end
        end
        
        axt(vv,cc) = subplot(length(VOIs),3,(vv-1)*3+cc);hold on;
        diff2plot = all_difft{vv,cc}; diff2plot(isnan(diff2plot))=0;
        imagesc(diff2plot);set(gca,'XTick',1:length(tpts_of_interest),'YTick',1:length(tpts_of_interest),'XTickLabel',tpts_of_interest*myTR,'YTickLabel',tpts_of_interest*myTR);
        colormap(rwb);
        %xlabel('T_2 (s)');
        %ylabel('T_1 (s)');
        if vv == 1
            title(condstr{cc});
        end
        
        if cc == 1
            ylabel(VOIs{vv});
        end
        set(gca,'XTickLabelRotation',90)
        
        axis square ij;
        
        
        
        
    end
 
    % draw square around relevant cells, black if FDR significant,
    % otherwise, (trend), gray square,
    tmpp = cell2mat({allp_difft{:,cc}});
    tmpp = tmpp(:);
    fdr_difft(cc) = fdr(tmpp(~isnan(tmpp)),0.05);
    % in case we want an FDR-thresh for each condition independently (ends
    % up being more conservative than across-condition concatenating before
    % FDR)
    
    
    fprintf('%s = FDR thresh = %0.03f\n',condstr{cc},fdr_difft(cc));
    
   
    
end

% FDR across all conditions together
tmpp = cell2mat({allp_difft{:}});
tmpp = tmpp(:);
fdr_difft_all = fdr(tmpp(~isnan(tmpp)),0.05);

for cc = 1:3
    for vv = 1:length(VOIs)
        axes(axt(vv,cc));
        for tt2 = 2:length(tpts_of_interest)
            
            for tt1 = 1:(tt2-1)

                if allp_difft{vv,cc}(tt2,tt1)<=fdr_difft_all
                    rectangle('Position',[(tt1-0.49) (tt2-0.49) 0.98 0.98],'EdgeColor',[0 0 0],'LineWidth',1,'Curvature',[0.1 0.1]);
                elseif allp_difft{vv,cc}(tt2,tt1)<=0.05
                    rectangle('Position',[(tt1-0.49) (tt2-0.49) 0.98 0.98],'EdgeColor',[0.5 0.5 00.5],'LineWidth',1,'Curvature',[0.1 0.1]);
                end  


            end
        end
    end
end
%match_clim(axt);
tmpclim = get(axt,'CLim');
cc = cell2mat({tmpclim{:}});

set(axt,'CLim',[-1 1]*max(abs(cc(:))),'XLim',[0.5 length(tpts_of_interest)+0.5],'YLim',[0.5 length(tpts_of_interest)+0.5])

if save_stats == 1
    fn_stats = sprintf('%swmDrop_stats/n%i_vectorMean_timeCourseDiffsOrig_%sthru%s_%s_%s.mat',root,length(u_subj),VOIs{1},VOIs{end},recon_str,datestr(now,30));
    fprintf('saving stats to %s...\n',fn_stats);
    save(fn_stats,'fdr_difft','fdr_difft_all','allp_difft','VOIs','n_resample_iter','rand_seed','delay_tpts');%,'mythresh_delay_compare','allp','allp_delay','allp_delay_compare','tpts_of_interest','VOIs','condstr','n_resample_iter','rand_seed','delay_tpts');
    %save(fn_stats,'fdr_difft','fdr_difft_all','allp_difft','allp_delay','allp_delay_compare','mythresh','mythresh_delay','mythresh_delay_compare','VOIs','n_resample_iter','rand_seed','delay_tpts');%,'mythresh_delay_compare','allp','allp_delay','allp_delay_compare','tpts_of_interest','VOIs','condstr','n_resample_iter','rand_seed','delay_tpts');
end



return


