function wmDrop_plotReconstructions_trnAvg_posCoreg1(subj,VOIs,tpts_of_interest, hex_size)
% plots WM reconstructions coregistered and sorted by relative position (Fig. 1B), for each TR
%
% Data presented in Fig. 4 of SES 2016 is generated using 51x51
% reconstructions, like follows:
%   wmDrop_plotReconstructions_posCoreg_thruTime1({'AI81','AI82','AI83','AP81','AP82','AP83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'},{'V1','V3A','IPS0','IPS2','sPCS','SuperWMDrop'});
%
% tpts_of_interest is in TR index (0 is TR closest to delay onset, 1 is 1
% TR after, etc). can try different delay period definitions, etc 
% use "hex_size" (along w/ channelRespAmp and computeReconstructions) to
% change # of channels in basis set (size and position are all yoked)


if nargin < 1    
    subj = {'AI81','AI82','AI83','AP81','AP82','AP83','AR81','AR82','AR83','AS81','AS82','AS83','AL81','AL82','AL83','BC81','BC82','BC83'};    
end

if nargin < 2

    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS','SuperWMDrop','iPCS','DLPFC','SMA'}; 

end

if nargin < 3
    tpts_of_interest = {[3 4],[7 8]};
end

if nargin < 4
    hex_size = 7;
end

root = load_root;



match_clim = 1;


u_subj = unique(cellfun(@(s) s(1:end-1),subj,'uniformoutput',0));


use_first_half_180 = 1;

condstr = {'R1','R2-neutral','R2-valid'};
targstr = {'PT','NPT'};

trials_per_superrun = 54;



nblank = length(subj) * length(VOIs) * 3 * trials_per_superrun * 12; % 3 superruns, 54 trials per superrun, 12 tpts

all_conds = nan(nblank,5);
all_tpts = nan(nblank,1);
all_subj = nan(nblank,1);
all_vois = nan(nblank,1);
all_180_half = nan(nblank,1);
all_coords{1} = nan(nblank,2);
all_coords{2} = nan(nblank,2);

% n_conds x n_VOIs
startidx = 1;
for ss = 1:length(subj)
    
    this_subj_id = find(strcmpi(u_subj,subj{ss}(1:end-1)));
    
    for vv = 1:length(VOIs)
        
        fn = sprintf('%swmDrop_recons/%s_%s_hex%i_trnAvg_pos_coreg1.mat',root,subj{ss},VOIs{vv},hex_size);
        fprintf('loading %s...\n',fn);
        load(fn);
        
        thisidx = startidx:(startidx+size(recons_vec{1},1)-1);
        
        if startidx == 1
            res = sqrt(size(recons_vec{1},2));
            
            
            [gridx,gridy] = meshgrid(linspace(-maxecc,maxecc,res),linspace(-maxecc,maxecc,res));
            gridx = reshape(gridx,numel(gridx),1);gridy = reshape(gridy,numel(gridy),1);
            
            
            all_recons{1} = nan(nblank,res^2);
            all_recons{2} = nan(nblank,res^2);
            
        end
        
        
        all_recons{1}(thisidx,:) = recons_vec{1};
        all_recons{2}(thisidx,:) = recons_vec{2};
        all_conds(thisidx,:) = conds;
        all_tpts(thisidx) = tpts;
        all_subj(thisidx) = this_subj_id*ones(size(conds,1),1);
        all_vois(thisidx) = vv*ones(size(conds,1),1);
        all_coords{1}(thisidx,:) = t_coord_cart{1};
        all_coords{2}(thisidx,:) = t_coord_cart{2};
        
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
        
        
        
        startidx = thisidx(end)+1;
        clear whichhalf
    end
end

valididx = 1:(startidx-1);
all_recons{1} = all_recons{1}(valididx,:);
all_recons{2} = all_recons{2}(valididx,:);
all_conds  = all_conds(valididx,:);
all_tpts = all_tpts(valididx);
all_subj = all_subj(valididx);
all_vois = all_vois(valididx);
all_180_half = all_180_half(valididx);


% like wmDelay Fig 3 - a figure for each condition/target, and a column
% within each for each separation condition (VOIs down each column)

su = unique(all_conds(:,3)); % 1, 2, 3 for 6 stim positions
for tt = 1:length(tpts_of_interest)
    ax{tt}= [];
    for ii = 2
        for cc = 1:3
            figure;
            for vv = 1:length(VOIs)
                
                for sepidx = 1:length(su)
                    
                    ax{tt}(end+1) = subplot(length(VOIs),length(su),sepidx+(vv-1)*length(su));
                    this_recon = nan(length(u_subj),size(all_recons{ii},2));
                    
                    % first average reconstructions for each subject, then average over
                    % subjects (so that different numbers of sessions don't contribute
                    % to differces in overall average) - this is how wmdelay plotting
                    % worked
                    for ss = 1:length(u_subj)
                        
                        if use_first_half_180==1
                            thisidx = all_subj==ss & all_conds(:,2)==cc & all_conds(:,3) == su(sepidx) & ismember(all_tpts,tpts_of_interest{tt}) & all_vois==vv & (isnan(all_180_half)|(all_180_half==1));
                        else
                            thisidx = all_subj==ss & all_conds(:,2)==cc & all_conds(:,3) == su(sepidx) & ismember(all_tpts,tpts_of_interest{tt}) & all_vois==vv;
                        end
                        this_recon(ss,:) = mean(all_recons{ii}(thisidx,:),1);
                        
                    end
                    
                    mm = mean(this_recon,1);
                    imagesc(reshape(mm,res,res));colormap parula;
                    
                    if vv == 1 && sepidx == ceil(length(su)/2)
                        title(sprintf('%s - %s',condstr{cc},targstr{ii}));
                    end
                    
                    if sepidx == 1
                        ylabel(VOIs{vv});
                    end
                    
                    axis square xy;
                    
                    clear this_recon;
                    
                end
            end
            if length(VOIs)==9
                set(gcf,'Position',[440    42   282   756]);
            end
        end
    end
    
    
end

if match_clim == 1
    tmpax = cell2mat(ax);
    thisax = tmpax(:);
    all_clim = cell2mat(get(thisax,'CLim'));
    set(thisax,'CLim',[min(all_clim(:,1)) max(all_clim(:,2))],'Box','on','XTick',[],'YTick',[]);
else
    
    for tt = 1:2
        thisax = ax{tt};
        all_clim = cell2mat(get(thisax,'CLim'));
        set(thisax,'CLim',[min(all_clim(:,1)) max(all_clim(:,2))],'Box','on','XTick',[],'YTick',[]);
        
    end
    
end

return
