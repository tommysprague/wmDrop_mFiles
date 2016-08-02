function wmDrop_plotReconstructions_posCoreg_thruTime1(subj,VOIs,tpts_of_interest,hex_size)
% plots WM reconstructions coregistered and sorted by relative position (Fig. 1B), for each TR
%
% Data presented in Fig. 3 of SES 2016 is generated using 101x101
% reconstructions, like follows:
%   wmDrop_plotReconstructions_posCoreg_thruTime1({'AI81','AI82','AI83','AP81','AP82','AP83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'},{'V1','V3A','IPS0','IPS2','sPCS','SuperWMDrop'});
%
% tpts_of_interest is in TR index (0 is TR closest to delay onset, 1 is 1 TR after, etc)
% use "hex_size" (along w/ channelRespAmp and computeReconstructions) to
% change # of channels in basis set (size and position are all yoked)

if nargin < 1
   
    subj = {'AI81','AI82','AI83','AP81','AP82','AP83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};

end

if nargin < 2
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS','SuperWMDrop','iPCS','DLPFC','SMA'}; % FEFnew
end

if nargin < 3
    tpts_of_interest = 0:9;
end

if nargin < 4
    hex_size = 7;
end


%custom_clim = [-0.25 0.4481];
custom_clim = [];

root = load_root;

myTR = 2.25;  % repetition time, in sec

%mycolors = [0 127 66; 56 9 124; 28 69 68]./255;
load([root 'wmDrop_colors.mat']);


u_subj = unique(cellfun(@(s) s(1:end-1),subj,'uniformoutput',0));

use_first_half_180 = 1;


condstr = {'R1','R2','R2d'};


trials_per_superrun = 54;

%tpts_of_interest = [2 3 4]; % to start with...

nblank = length(subj) * length(VOIs) * 3 * 54 * length(tpts_of_interest); % 3 superruns, 54 trials per superrun, 10 tpts
%all_recons = nan(nblank,res^2);
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
        
        fn = sprintf('%swmDrop_recons/%s_%s_hex%i_trnAvg_pos_coreg1.mat',root,subj{ss},VOIs{vv},hex_size);
        load(fn);
        
        thisidx = startidx:(startidx+size(recons_vec{1},1)-1);
        
        if startidx == 1
            res = sqrt(size(recons_vec{1},2));
            
            
            [gridx,gridy] = meshgrid(linspace(-maxecc,maxecc,res),linspace(-maxecc,maxecc,res));
            gridx = reshape(gridx,numel(gridx),1);gridy = reshape(gridy,numel(gridy),1);

            
            all_recons = nan(nblank,res^2);
        end
        
        all_recons(thisidx,:) = recons_vec{2};
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
                tmp180 = find(tmpconds(:,3)==3&tmpconds(:,2)==rr);
                whichhalf(sridx(tmp180(1:(length(tmp180)/2)))) = 1;
                whichhalf(sridx(tmp180( (1+(length(tmp180)/2)):length(tmp180)))) = 2;
                %all_sridx = thisidx(sridx); % this is where to put the values in
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


ax = [];

% target separation condition plotted in Fig. 3 is 2 - try 1 or 3 as well
sep_of_interest = 2;

% loop over conditions (1 = R1, 2 = R2-neutral, 3 = R2-valid)
for cc = 1:3    
    figure;
    for tt = 1:length(tpts_of_interest)
        for vv = 1:length(VOIs)
            ax(end+1) = subplot(length(VOIs),length(tpts_of_interest),tt+(vv-1)*length(tpts_of_interest));
            hold on;
            this_recon = nan(length(u_subj),size(all_recons,2));
            
            for ss = 1:length(u_subj)
                
                
                if use_first_half_180 == 1
                    thisidx = all_subj==ss & all_conds(:,2)==cc & all_conds(:,3)==sep_of_interest & ismember(all_tpts,tpts_of_interest(tt)) & all_vois==vv & (isnan(all_180_half)|(all_180_half==1));
                else
                    thisidx = all_subj==ss & all_conds(:,2)==cc & all_conds(:,3)==sep_of_interest & ismember(all_tpts,tpts_of_interest(tt)) & all_vois==vv;
                end
                
                this_recon(ss,:) = mean(all_recons(thisidx,:),1);
                
            end
            mm = mean(this_recon,1);
            imagesc(gridx,gridy,reshape(mm,res,res));colormap parula;
            
            if vv == 1
                % to do: in color of cc...
                title(sprintf('%.02f s',tpts_of_interest(tt)*myTR));
            end
            
            if tt == 1
                ylabel(VOIs{vv});
            end
            
            if vv == length(VOIs) && tt == length(tpts_of_interest)   % draw a 3 dva scalebar in bottom left
                plot([-5.5 -2.5],[-5.0 -5.0],'k-','LineWidth',1);
            end
            
            axis square xy tight;
            hold off;
            clear this_recon;
        end
    end


    set(gcf,'Position',[440   460   933   338],'Name',['Fig. 3 - ' condstr{cc}],'NumberTitle','off');
    
end
set(ax,'XTick',[],'YTick',[],'Box','on');
if ~isempty(custom_clim)
    set(ax,'CLim',custom_clim);
else
    all_clim = cell2mat(get(ax,'CLim'));
    set(ax,'CLim',[min(all_clim(:,1)) max(all_clim(:,2))],'XTick',[],'YTick',[]);
end
return
