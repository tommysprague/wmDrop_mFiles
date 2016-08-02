function wmDrop_plotFits_trnAvg_exactCoreg_byRecallError1_resamp(n_subj,VOIs,tpts_of_interest,hex_size)
%
% plots fits to each half of data, split by median recall error, separately
% loads files created with wmDrop_fitReconstructions_exactCoreg_byRecallError1_resample.m
% (ALLn3_V1_R1_tpts34_ms1_hex7_Iter1000_fitExactCoreg1_resample.mat)
%
% For Sprague, Ester & Serences, 2016, Fig. 8, plot only SuperWMDrop; for
% Fig. S8, plot V1-sPCS. Stats computed that way as well.

save_stats = 1;
if nargin < 1
    n_subj = 6;
end

if nargin < 2
    VOIs = {'V1','V2','V3','V3A','V4','IPS0','IPS1','IPS2','IPS3','sPCS'};
    
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
ms_str = {'Low error','High error'};

% variables to load
nblank = length(VOIs) * length(r_cond) * length(tpts_of_interest) * 2;

all_cond = nan(nblank,1);
all_tpts = nan(nblank,1);
all_VOIs = nan(nblank,1);

all_medianSplit = nan(nblank,1);
all_exflag = nan(nblank,n_iter);


startidx = 1;

for ms = 1:2
    for rr = 1:length(r_cond)
        for vv = 1:length(VOIs)
            for tt = 1:length(tpts_of_interest)
                
                myfn = sprintf('%swmDrop_fits/ALLn%i_%s_R%i_tpts%s_ms%i_hex%i_Iter%i_fitExactCoreg_byRecallError1_CR.mat',root,n_subj,VOIs{vv},r_cond(rr),tptstr{tt},ms,hex_size,n_iter);
                fprintf('loading %s...\n',myfn);
                fits = load(myfn); clear myfn;
                
                if vv == 1 && tt == 1 && rr == 1 && ms == 1
                    all_recon = nan(nblank,size(fits.m_recon_vec,2));
                    all_bf = nan(nblank,size(fits.bf_fine,2),n_iter);
                end
                
                thisidx = startidx;
                
                % save labels
                all_cond(thisidx) = rr;
                all_tpts(thisidx) = tt;
                
                all_VOIs(thisidx) = vv;
                all_medianSplit(thisidx) = ms;
                
                % save fit params/recons
                all_recon(thisidx,:) = mean(fits.m_recon_vec,1);
                all_bf(thisidx,:,:) = fits.bf_fine'; % FINE
                
                all_exflag(thisidx,:) = fits.ex_flag_fine;
                
                startidx = thisidx(end)+1;
                
                clear fits thisidx;
            end
        end
    end
end

% put size paramter into fwhm units
all_bf(:,3,:) = rad2fwhm(all_bf(:,3,:));


axi = [];
th = linspace(0,2*pi,1001);

% plot mean reconstructions w/ mean fit (as dotted circle, like
% before)
% 4 figures (for 2 delays, 2 median splits), 3 rows, n_vois columns
for tt = 1:length(tpts_of_interest)
    
    for ms = 1:2
        
        figure;
        for vv = 1:length(VOIs)
            
            
            
            for cc = 1:length(condstr)
                
                thisidx = all_VOIs == vv & all_cond==cc & all_tpts==tt & all_medianSplit == ms;
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
            
            
            set(gcf,'Position',[440   365   463   433],'NumberTitle','off','Name',sprintf('%s - %s',delay_str{tt},ms_str{ms}));
            
            
            
        end
    end
end

all_clim = cell2mat(get(axi,'CLim'));
set(axi,'CLim',[min(all_clim(:,1)) max(all_clim(:,2))]);


% plot one subplot for size, amplitude, baseline for each delay period
% in each subplot, for each ROI break up low and high recall error

relx = linspace(-.15,.15,length(ms_str));

ax = nan(length(params_of_interest),length(tpts_of_interest)*length(condstr));
for cc = 1:length(condstr)
    figure;
    for pp = 1:length(params_of_interest)
        
        for tt = 1:length(tpts_of_interest)
            ax(pp,(tt-1)*length(condstr)+cc) = subplot(length(params_of_interest),length(tpts_of_interest),(pp-1)*length(tpts_of_interest)+tt); hold on;
            for vv = 1:length(VOIs)
                
                
                for ms = 1:2
                    
                    thisidx = all_VOIs == vv & all_cond == cc & all_tpts==tt & all_medianSplit == ms;
                    
                    mm = mean(all_bf(thisidx,params_of_interest(pp),:),3);
                    ci = prctile(squeeze(all_bf(thisidx,params_of_interest(pp),:)),[2.5 97.5]);
                    plot(vv+relx(ms)*[1 1],ci,'-','linewidth',1,'color',mycolors(cc,:));clear mm ci;
                    
                    if ms == 1
                        plot(vv+relx(ms),mean(all_bf(thisidx,params_of_interest(pp),:),3),'o','color',mycolors(cc,:),'MarkerFaceColor',mycolors(cc,:),'LineWidth',1,'MarkerSize',7);
                    else
                        plot(vv+relx(ms),mean(all_bf(thisidx,params_of_interest(pp),:),3),'o','color',mycolors(cc,:),'MarkerFaceColor',[1 1 1],'LineWidth',1,'MarkerSize',7);
                    end
                    
                    
                end
                
            end
            if tt == 1
                ylabel(param_str{pp});
            end
            if pp == 1
                title(delay_str{tt});
            end
        end
        
        
    end
end

for pp = 1:length(params_of_interest)
    tmpylim = cell2mat(get(ax(pp,:),'YLim'));
    set(ax(pp,:),'YLim',[min(tmpylim(:,1)) max(tmpylim(:,2))]);
    
end
set(ax,'XLim',[0 length(VOIs)+1]);
set(ax,'XTick',1:length(VOIs),'XTickLabel',VOIs,'TickLength',[0.025 0.025],'box','off');


%% compute p-values between median splits

allp_medsplit = cell(length(params_of_interest),1);
fdr_medsplit = nan(length(params_of_interest),1);

for pp = 1:length(params_of_interest)
    allp_medsplit{pp} = nan(length(condstr),length(VOIs),length(tpts_of_interest));
    for cc = 1:length(condstr)
        
        % for now, doing all these on same figure - can break into two figs if
        % necessary
        for tt = 1:length(tpts_of_interest)
            for vv = 1:length(VOIs)
                fprintf('Cond: %s\tParam: %s\t Delay: %i\tROI: %s\n',condstr{cc},param_str{pp},tt,VOIs{vv});
                
                thisidx1 = all_VOIs == vv & all_cond == cc & all_tpts==tt & all_medianSplit == 1;
                thisidx2 = all_VOIs == vv & all_cond == cc & all_tpts==tt & all_medianSplit == 2;
                
                diffs = squeeze(all_bf(thisidx1,params_of_interest(pp),:)) - squeeze(all_bf(thisidx2,params_of_interest(pp),:));
                thisp = 2*min(mean(diffs<0),mean(diffs>0));
                allp_medsplit{pp}(cc,vv,tt) = thisp;
                fprintf('p = %0.03f\n\n',thisp);
            end
        end
    end
    
    fdr_medsplit(pp) = fdr(allp_medsplit{pp}(:),0.05);
    
end

if length(VOIs) > 1
    voi_str = [VOIs{1} 'to' VOIs{end}];
else
    voi_str = VOIs{1};
end

if save_stats == 1
    fn2s = sprintf('%swmDrop_stats/n%i_plotFits_exactCoreg_byRecallError_%s_%s.mat',root,n_subj,voi_str,datestr(now,30));
    
    fprintf('saving to %s...\n',fn2s);
    save(fn2s,'allp_medsplit','fdr_medsplit','n_subj','VOIs','params_of_interest','condstr');
end
return