function wmDrop_plotFits_trnAvg_exactCoreg1_resample(n_subj,VOIs,tpts_of_interest,hex_size)
%
% WMDROP_PLOTFITS_TRNAVG_EXACTCOREG1_RESAMPLE 
% plots fits and resampled reconstructions (with fits) as in Fig. 7 of
% Sprague, Ester & Serences, 2016. Note that we only include n_subj here,
% as all fits are saved after resampling across subjects. this all may
% require some modification if only particular subj are to be plotted
%
% for printed stats (bottom), must do 10 individual ROIs and Super-ROI
% separately. FDR correction is based on group of ROIs run; because
% SuperROI is not independent of any others, cannot be subjected to same
% corrections. 

save_stats = 1;

if nargin < 1
    n_subj = 6;
end

if nargin < 2 
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS','SuperWMDrop'};

end


r_cond = [1 2 3];

if nargin < 3
    tpts_of_interest = {[3 4],[7 8]};
end

if nargin < 4
    hex_size = 7;
end

if ~iscell(tpts_of_interest)
    tpts_of_interest{1} = tpts_of_interest;
end

for tt = 1:length(tpts_of_interest)
    tptstr{tt} = '';
    for ii = 1:length(tpts_of_interest{tt})
        tptstr{tt} = [tptstr{tt} num2str(tpts_of_interest{tt}(ii))];
    end
end

root = load_root;

n_iter = 1000;

params_of_interest = [3 4 5];
param_str = {'Size','Amplitude','Baseline'};

delay_str = {'First delay','Second delay'};


ecc_to_align = 3.5; % align to this point


load(fullfile(root,'wmDrop_colors.mat'));

condstr = {'R1','R2-neutral','R2-valid'};
load(sprintf('%swmDrop_subjSym.mat',root)); 

% variables to load
nblank = length(VOIs) * length(r_cond) * length(tpts_of_interest);
all_cond = nan(nblank,1);
all_tpts = nan(nblank,1);
all_VOIs = nan(nblank,1);   
all_subj = nan(nblank,1);
all_exflag = nan(nblank,n_iter);


startidx = 1;

%for ss = 1:length(subj)
for rr = 1:length(r_cond)
    for vv = 1:length(VOIs)
        for tt = 1:length(tpts_of_interest)
            
            myfn = sprintf('%swmDrop_fits/ALLn%i_%s_R%i_tpts%s_hex%i_Iter%i_fitExactCoreg1_CR.mat',root,n_subj,VOIs{vv},r_cond(rr),tptstr{tt},hex_size,n_iter);
            
            fprintf('loading %s...\n',myfn);
            fits = load(myfn); clear myfn;            
            
            if vv == 1 && tt == 1 && rr == 1
                all_recon = nan(nblank,size(fits.m_recon_vec,2));
                all_bf = nan(nblank,size(fits.bf_fine,2),n_iter);
                all_bfg = nan(nblank,size(fits.bf_grid,2),n_iter);

            end
            
            thisidx = startidx;%startidx:(startidx+size(fits.cond,1)-1);
            
            % save labels
            all_cond(thisidx) = rr;
            all_tpts(thisidx) = tt;
            all_VOIs(thisidx) = vv;
            
            % save fit params/recons
            all_recon(thisidx,:) = mean(fits.m_recon_vec,1);
            all_bf(thisidx,:,:) = fits.bf_fine'; % FINE
            all_bfg(thisidx,:,:) = fits.bf_grid';
            all_exflag(thisidx,:) = fits.ex_flag_fine;
            
            startidx = thisidx(end)+1;
            
            clear fits thisidx;
        end
    end
end

% put size paramter into fwhm units
all_bf(:,3,:) = rad2fwhm(all_bf(:,3,:));


axi = [];
th = linspace(0,2*pi,1001);

% plot mean reconstructions w/ mean fit (as dotted circle, like
% before)
% 2 figures (for 2 delays), 3 rows, n_vois columns 
for tt = 1:length(tpts_of_interest)
    
    figure;  
    for vv = 1:length(VOIs)
        
        
        
        for cc = 1:length(condstr)
            
            thisidx = all_VOIs == vv & all_cond==cc & all_tpts==tt;
            thisrecon = mean(all_recon(thisidx,:),1);
            
            thisax=subplot(length(r_cond),length(VOIs),(cc-1)*length(VOIs)+vv); hold on;
            imagesc([-6 6],[-6 6],reshape(thisrecon,sqrt(length(thisrecon)),sqrt(length(thisrecon))));
            
            axis xy square;
            
            plot(ecc_to_align,0,'o','MarkerFaceColor',[0.8500    0.3250    0.0980],'Color',[0.8500    0.3250    0.0980],'MarkerSize',5);
            
            xx = mean(all_bf(thisidx,3,:),3)*0.5*cos(th)+mean(all_bf(thisidx,1,:),3);
            yy = mean(all_bf(thisidx,3,:),3)*0.5*sin(th)+mean(all_bf(thisidx,2,:),3);
            plot(xx,yy,'w--','LineWidth',1.5);
            plot(mean(all_bf(thisidx,1,:),3),mean(all_bf(thisidx,2,:),3),'w+','MarkerSize',3,'LineWidth',2);
            
            
            if cc == 1
                title(VOIs{vv});
            end
            
            if vv == 1
                ylabel(condstr{cc});
            end
            
            set(thisax,'XTick',[-ecc_to_align 0 ecc_to_align],'YTick',[-ecc_to_align 0 ecc_to_align],'TickLength',[0.025 0.025],'XTickLabel',[],'YTickLabel',[],'XLim',[-6 6],'YLim',[-6 6],'box','on');
            
            
            
            axi(end+1) = thisax; clear thisax thisidx thisrecon;
        end
        
        
        set(gcf,'Position',[440   365   463   433],'NumberTitle','off','Name',sprintf('%s',delay_str{tt}));
        

        
    end
end

all_clim = cell2mat(get(axi,'CLim'));
set(axi,'CLim',[min(all_clim(:,1)) max(all_clim(:,2))]);


% plot one subplot for size, amplitude, baseline for each delay period
% in each subplot, have a group of dots for each ROI, and each subject and
% each condition as separate columns of dots (like Fig 4D of CB - but
% include individual subject datapoints)

relx = linspace(-.2,.2,length(condstr));

ax = nan(length(params_of_interest),length(tpts_of_interest));

figure;
for pp = 1:length(params_of_interest)

    
    for tt = 1:length(tpts_of_interest)
        ax(pp,tt) = subplot(length(params_of_interest),length(tpts_of_interest),(pp-1)*length(tpts_of_interest)+tt); hold on;
        for vv = 1:length(VOIs)
            
            
            for cc = 1:length(condstr)

                thisidx = all_VOIs == vv & all_cond == cc & all_tpts==tt;
                plot(vv+relx(cc),mean(all_bf(thisidx,params_of_interest(pp),:),3),'o','color',mycolors(cc,:),'MarkerFaceColor',mycolors(cc,:),'LineWidth',1,'MarkerSize',7);
                mm = mean(all_bf(thisidx,params_of_interest(pp),:),3);
                ci = prctile(squeeze(all_bf(thisidx,params_of_interest(pp),:)),[2.5 97.5]);%std(all_bf(thisidx,params_of_interest(pp)))/sqrt(length(subj));
                plot(vv+relx(cc)*[1 1],ci,'-','linewidth',1,'color',mycolors(cc,:));clear mm ci;
            end
            
        end
        if tt == 1
            ylabel(param_str{pp});
        end
        if pp == 1
            title(delay_str{tt});
        end
    end
    tmpylim = cell2mat(get(ax(pp,:),'YLim'));
    set(ax(pp,:),'YLim',[min(tmpylim(:,1)) max(tmpylim(:,2))]);
end
set(ax,'XLim',[0 length(VOIs)+1]);
set(ax,'XTick',1:length(VOIs),'XTickLabel',VOIs,'TickLength',[0.02 0.02],'box','off');


%% print stats

cond_comparisons = [1 2; 2 3; 1 3];  % compare col 1 & col 2 

allp_tab = cell(length(params_of_interest),1);

for pp = 1:length(params_of_interest)
    allp = [];
    allp_tab{pp} = nan(length(tpts_of_interest)*length(VOIs),3);
    pidx = 1;
    for tt = 1:length(tpts_of_interest)
        for vv = 1:length(VOIs)
            fprintf('%s\tDelay %i\t%s\n',param_str{pp},tt,VOIs{vv});
            for cc = 1:3
                thisidx1 = all_VOIs == vv & all_cond == cond_comparisons(cc,1) & all_tpts==tt;
                thisidx2 = all_VOIs == vv & all_cond == cond_comparisons(cc,2) & all_tpts==tt;
                mydiff = squeeze(all_bf(thisidx1,params_of_interest(pp),:)) - squeeze(all_bf(thisidx2,params_of_interest(pp),:));
                myp = min(mean(mydiff<0),mean(mydiff>0))*2;
                allp = [allp;myp];
                fprintf('%s\tvs\t%s: p = %0.03f\n',condstr{cond_comparisons(cc,1)},condstr{cond_comparisons(cc,2)},myp);
                
                allp_tab{pp}(pidx,cc) = myp;
                
                clear myp
            end
            pidx = pidx+1;
            
        end
    end
    fdr_thresh(pp) = fdr(allp,0.05); clear allp;
    fprintf('%s - FDR threshold: p = %0.03f\n',param_str{pp},fdr_thresh(pp));
end

if save_stats == 1
    if length(VOIs)==1
        fn2s = sprintf('%swmDrop_stats/n%i_fits_exactCoreg1_CR_%iIter_%s_%s.mat',root,n_subj,n_iter,VOIs{1},datestr(now,30));
    else
        fn2s = sprintf('%swmDrop_stats/n%i_fits_exactCoreg1_CR_%iIter_%iVOIs_%s.mat',root,n_subj,n_iter,length(VOIs),datestr(now,30));
    end
    fprintf('saving stats to %s\n',fn2s);
    save(fn2s,'fdr_thresh','allp_tab','n_subj','VOIs','tpts_of_interest','cond_comparisons');
end
return