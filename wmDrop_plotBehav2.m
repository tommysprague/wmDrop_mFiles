function wmDrop_plotbehav2(subj)

% TCS 11/19/2015

close all;

% this seed used for publication
rng(234343);

root = load_root;%'/usr/local/serenceslab/tommy/wmDrop/';

load([root 'wmDrop_colors.mat']);
load([root 'wmDrop_subjSym.mat'])


if nargin < 1
    subj = {'AI81','AI82','AI83','AP81','AP82','AP83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};
end

u_subj = unique(cellfun(@(s) s(1:end-1),subj,'uniformoutput',0));


n_scanner_trials = 3*54; % maximum number of trials per session

nblank = n_scanner_trials*length(subj);

all_conds = nan(nblank,5);
all_resp  = nan(nblank,1);
all_respErr = nan(nblank,1);
all_TRcoord = nan(nblank,2);
all_TNcoord = nan(nblank,2);

all_respCoord = nan(nblank,1);


all_subj = nan(nblank,1);

startidx = 1;
for ss = 1:length(subj)
    
    this_subj_id = find(strcmpi(u_subj,subj{ss}(1:end-1)));
    
    fn = sprintf('%swmDrop_behav/%s_wmDrop1dScanner_behav.mat',root,subj{ss});
    fprintf('loading %s...\n',fn);
    thisfile = load(fn);
    
    thisidx = startidx:(startidx+size(thisfile.all_conds,1)-1);
    
    all_conds(thisidx,:) = thisfile.all_conds;
    all_TRcoord(thisidx,:) = thisfile.all_TRcoord;
    all_TNcoord(thisidx,:) = thisfile.all_TNcoord;
    
    all_resp(thisidx) = thisfile.all_resp;
    all_respErr(thisidx) = thisfile.all_respErr;
    
    all_respCoord(thisidx) = thisfile.all_respCoord;
    
    all_subj(thisidx) = this_subj_id;
    
    
    startidx = thisidx(end)+1;
    clear thisfile
end



allidx = 1:(startidx-1);
all_conds = all_conds(allidx,:);
all_TRcoord = all_TRcoord(allidx,:);
all_TNcoord = all_TNcoord(allidx,:);
all_resp = all_resp(allidx);
all_respErr = all_respErr(allidx);
all_subj = all_subj(allidx);

all_respCoord = all_respCoord(allidx);


include_trial = ones(startidx-1,1)==1; % for now, but later we can drop trials based on artifact rejection; behav performance, etc


%% mean recall error for R1trn, R1tst, and R2d
% plot mean as bar graph, then individual participants in front

conds_of_interest = [1;2;3]; % 1st col & 2nd col
cond_str = {'R1','R2','R2-valid'};
% all_meanRecallErr
all_mre = nan(length(u_subj),size(conds_of_interest,1));


for cc = 1:size(conds_of_interest,1)
    %subplot(1,size(conds_of_interest,1),cc);
    for ss = 1:length(u_subj)
        thisidx = all_conds(:,2)==conds_of_interest(cc) & all_subj==ss & include_trial==1;
        all_mre(ss,cc) = mean(all_respErr(thisidx));
        %all_mre(ss,cc) = median(all_respErr(thisidx));
        clear thisidx;
    end
end



%% alternative means of visualizing (bars instead of CIs)

figure; hold on;
% bars;
for cc = 1:size(conds_of_interest,1)
    bar(cc,mean(all_mre(:,cc),1),'FaceColor',mycolors(conds_of_interest(cc),:));
end

% subj lines
for ss = 1:length(u_subj)
    for cc = 1:(size(conds_of_interest,1)-1)
        xidx = [cc cc+1];
        plot(xidx,all_mre(ss,xidx),'-','Color',[0.5 0.5 0.5]);clear xidx;
    end
    for cc = 1:size(conds_of_interest,1)
        plot(cc,all_mre(ss,cc),subjsym.(u_subj{ss}),'Color',mycolors(conds_of_interest(cc),:),'MarkerFaceColor',[1 1 1],'MarkerSize',8);
    end
end

set(gca,'XTick',[1 2 3],'XTickLabel',{'R1','R2','R2-valid'});
ylabel('Recall error (\circ)');
xlabel('Condition');
title('Mean recall error');
hold off;








%% stats

n_resample_iter = 1000;

% all_meanRecallErr
all_mre_resample = nan(size(conds_of_interest,1),n_resample_iter);
resample_idx = cell(size(conds_of_interest,1),1);%,n_resample_iter);



for cc = 1:size(conds_of_interest,1)
    
    thisidx = all_conds(:,2)==conds_of_interest(cc) & include_trial==1 ;
    thisdata = all_respErr(thisidx);
    
    resample_idx{cc} = nan(n_resample_iter,sum(thisidx));
    
    for ii = 1:n_resample_iter
        
    
        resample_idx{cc}(ii,:) = randsample(length(thisdata),length(thisdata),true);
        
        
        
        all_mre_resample(cc,ii) = mean(thisdata(resample_idx{cc}(ii,:)));
        
    end
    clear thisidx thisdata;
end

cond_compare = [1 2; 2 3; 1 3];
for cc = 1:size(cond_compare,1)
    thisp = 2*min( mean((all_mre_resample(cond_compare(cc,1),:)-all_mre_resample(cond_compare(cc,2),:)) < 0),mean((all_mre_resample(cond_compare(cc,1),:)-all_mre_resample(cond_compare(cc,2),:)) > 0));
    fprintf('R%i vs R%i: p = %0.03f\n',cond_compare(cc,1),cond_compare(cc,2),thisp);
end


%% Fig1c for Neuron - now has CIs, means

box_width = 0.35; % distance left, right of box center

% copied from scanner-only above
figure; hold on;
% boxes:
for cc = 1:size(conds_of_interest,1)
    thisci = prctile(all_mre_resample(cc,:),[2.5 97.5]);
    % draw a white box with colored borders
    r = rectangle('Position',[cc-box_width thisci(1) 2*box_width (thisci(2)-thisci(1))]);
    set(r,'FaceColor',[1 1 1],'EdgeColor',mycolors(conds_of_interest(cc),:),'LineWidth',1);
    
    % then a line at its mean
    line([cc-box_width cc+box_width].',mean(all_mre_resample(cc,:)) * [1 1].','LineWidth',1,'Color',mycolors(conds_of_interest(cc),:));
    %bar(cc,mean(all_mre(:,cc),1),'FaceColor',mycolors(conds_of_interest(cc),:));
end

% subj lines
for ss = 1:length(u_subj)
    for cc = 1:(size(conds_of_interest,1)-1)
        xidx = [cc cc+1];
        plot(xidx,all_mre(ss,xidx),'-','Color',[0.5 0.5 0.5]);clear xidx;
    end
    for cc = 1:size(conds_of_interest,1)
        plot(cc,all_mre(ss,cc),subjsym.(u_subj{ss}),'Color',mycolors(conds_of_interest(cc),:),'MarkerFaceColor',[1 1 1],'MarkerSize',8);
    end
end

set(gca,'XTick',[1 2 3],'XTickLabel',{'R1','R2','R2\nvalid'});
ylabel('Recall error (\circ)');
xlabel('Condition');
title('Mean recall error');
hold off;


%% Error histograms (Fig. S1 of Neuron)
% one column each subj, one row each condition


mybins = linspace(0,max(all_respErr),21);

meanindicator = [];

figure;ax = [];
for ss = 1:length(u_subj)
    for cc = 1:size(conds_of_interest,1)
        
        ax(end+1) = subplot(size(conds_of_interest,1),length(u_subj),(cc-1)*length(u_subj)+ss);
        hold on;
        
        
        
        thisidx = all_conds(:,2)==conds_of_interest(cc) & all_subj==ss & include_trial==1 ;
        
        
        this_hist = hist(all_respErr(thisidx),mybins);
        this_hist = this_hist/sum(this_hist);
        bar(mybins,this_hist,'BarWidth',1.75*(mybins(2)-mybins(1)),'FaceColor',mycolors(conds_of_interest(cc),:),'EdgeColor',mycolors(conds_of_interest(cc),:));
        
        %meanindicator(end+1) = plot(mean(all_respErr(thisidx)),0.10,'kv','MarkerFaceColor',[0 0 0]);
        meanindicator(end+1) = plot([1 1]*mean(all_respErr(thisidx)),[0 0.05],'k-');
        
        if cc == 1
            title(sprintf('%s',u_subj{ss}));
        end
        
        if ss == 1
            ylabel(cond_str{cc});
        else
            set(gca,'YTickLabel',[]);
        end
        
        if ss == 1 && cc == size(conds_of_interest,1)
            xlabel('Recall error (\circ)');
        else
            set(gca,'XTickLabel',[]);
        end
        
        hold off;
        clear this_hist;
    end
end

%match_ylim(ax);

set(ax,'XLim',[-1*(mybins(2)-mybins(1)) max(mybins)+(mybins(2)-mybins(1))],'YTick',[0:.5:1],'YLim',[0 0.5],'TickLength',[0.0175 0.0175]);
tmp = get(gca,'YLim');
set(meanindicator,'YData',[0 tmp(2)]);


return
