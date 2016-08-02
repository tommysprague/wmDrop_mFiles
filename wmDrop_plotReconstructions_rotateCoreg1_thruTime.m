function wmDrop_plotReconstructions_rotateCoreg1_thruTime(subj,VOIs,tpts_of_interest, hex_size)

% plots 1d reconstructions through time (Fig. 5)
%
% T Sprague (tommy.sprague@gmail.com)

if nargin < 1
    subj = {'AP81','AP82','AP83','AI81','AI82','AI83','AL81','AL82','AL83','AR81','AR82','AR83','AS81','AS82','AS83','BC81','BC82','BC83'};
    

end

if nargin < 2
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS','SuperWMDrop'};

end

if nargin < 3
    tpts_of_interest = 0:9;
end

if nargin < 4
    hex_size = 7;
end

custom_clim = [];

root = load_root;



myTR = 2.25;  % repetition time, in sec

load(fullfile(root, 'wmDrop_colors.mat'));



u_subj = unique(cellfun(@(s) s(1:end-1),subj,'uniformoutput',0));


use_first_half_180 = 1;

condstr = {'R1','R2-neutral','R2-valid'};
recon_str = 'trnAvg1'; 


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
        
        if ss == 1 && vv == 1
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

valididx = 1:(startidx-1);
all_recons = all_recons(valididx,:);
all_conds  = all_conds(valididx,:);
all_tpts = all_tpts(valididx);
all_subj = all_subj(valididx);
all_vois = all_vois(valididx);
all_180_half = all_180_half(valididx);



%% plot timecourse image (tpts rows x rad_pix cols) - like Eddie's figures


gridt = reshape(gridt,res_r,res_t);
gridr = reshape(gridr,res_r,res_t);

myth = gridt(1,:); % for x axis of reconstruciton plots
myth = 360*myth/(2*pi);
myr = gridr(:,1).';
which_conds = [1 2 3];

ax = [];
%
figure;

for vv = 1:length(VOIs)
    for cc = 1:length(which_conds)
        
        thisax = subplot(length(VOIs),length(which_conds),cc+(vv-1)*length(which_conds));
        ax(end+1)=thisax;
        hold on;
        
        this_recon = nan(length(tpts_of_interest),length(myth),length(u_subj));
        
        
        % now fill up this_recon
        for tt = 1:length(tpts_of_interest)
            for ss = 1:length(u_subj)
                if use_first_half_180==1
                    thisidx = all_subj==ss & all_conds(:,2)==cc & all_tpts==tpts_of_interest(tt) & all_vois==vv & (isnan(all_180_half)|(all_180_half==1));
                else
                    thisidx = all_subj==ss & all_conds(:,2)==cc & all_tpts==tpts_of_interest(tt) & all_vois==vv;
                end
                
                tmprecon = mean(all_recons(thisidx,:));
                tmprecon = reshape(tmprecon,res_r,res_t);
                this_recon(tt,:,ss) = mean(tmprecon,1);
                clear tmprecon;
                
            end
        end
        
        
        mm = mean(this_recon,3);
        
        imagesc(tpts_of_interest*myTR,myth,mm.');colormap parula;
        
        axis tight;
        set(thisax,'XTick',0:4:(myTR*tpts_of_interest(end)),'YTick',[-90 0 90],'TickLength',[0.025 0.025],'Box','on','YDir','reverse');%'XLim',[min(myth) max(myth)],...
        
        if vv == 1
            title(condstr{cc});
        end
        
        if vv~=length(VOIs)
            set(thisax,'XTickLabel',[]);
        end
        
        
        if cc == 1
            ylabel(VOIs{vv});
            
        end
        
        if cc~=1
            set(thisax,'YTickLabel',[]);
        end
       
        
        hold off;
        clear this_recon;
        
    end
    if length(u_subj)==2
        set(gcf,'Position',[440   635   956   163]);
    end
end
if ~isempty(custom_clim)
    set(ax,'CLim',custom_clim);
else
    all_clim = cell2mat(get(ax,'CLim'));
    set(ax,'CLim',[min(all_clim(:,1)) max(all_clim(:,2))]);
end

return

