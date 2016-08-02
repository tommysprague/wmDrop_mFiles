function wmDrop_plotReconstructions_rotateCoreg1_resample(n_subj,VOIs,tpts_of_interest,hex_size)
%
% plots 1-d resampled reconstructions w/ CIs
% 




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

lim_val = 0.02;


delay_str = {'First delay','Second delay'};


recon_str = 'trnAvg1'; 


load(fullfile(root,'wmDrop_colors.mat'));



condstr = {'R1','R2-neutral','R2-valid'};

% variables to load
nblank = length(VOIs) * length(r_cond) * length(tpts_of_interest);
all_cond = nan(nblank,1);
all_tpts = nan(nblank,1);
all_VOIs = nan(nblank,1);   



startidx = 1;

%for ss = 1:length(subj)
for rr = 1:length(r_cond)
    for vv = 1:length(VOIs)
        for tt = 1:length(tpts_of_interest)
            
         
            myfn = sprintf('%swmDrop_recons/ALLn%i_%s_R%i_tpts%s_hex%i_Iter%i_%s_resampleRotateCoreg1.mat',root,n_subj,VOIs{vv},r_cond(rr),tptstr{tt},hex_size,n_iter,recon_str);
            
            fprintf('loading %s...\n',myfn);
            reconstructions = load(myfn); clear myfn;            
            
            if vv == 1 && tt == 1 && rr == 1
                all_recon = nan(nblank,size(reconstructions.m_recon_vec,2));
                all_ci = nan(nblank,size(reconstructions.m_recon_vec,2),2);
         
                myth = reconstructions.myt.*360./(2*pi);
            end
            
            thisidx = startidx;
            
            % save labels
            all_cond(thisidx) = rr;
            all_tpts(thisidx) = tt;
            
            all_VOIs(thisidx) = vv;
            
            % save fit params/recons
            all_recon(thisidx,:) = mean(reconstructions.m_recon_vec,1);
            all_ci(thisidx,:,1) = prctile(reconstructions.m_recon_vec,[2.5],1);
            all_ci(thisidx,:,2) = prctile(reconstructions.m_recon_vec,[97.5],1);
            
            startidx = thisidx(end)+1;
            
            clear fits thisidx;
        end
    end
end




%% plot resampled means of 1-d reconstructions w/ CIs

axi = [];
%th = linspace(-pi,pi,1001);
% column-wise! each row a VOI, each column a condition
% 
% 2 figures (for 2 delays) 
for tt = 1:length(tpts_of_interest)
    
    figure;  
    for vv = 1:length(VOIs)
        
        
        
        for cc = 1:length(condstr)
            
            thisidx = all_VOIs == vv & all_cond==cc & all_tpts==tt;
            thisrecon = mean(all_recon(thisidx,:),1);
            thisci = squeeze(all_ci(thisidx,:,:)).';
            
            
            thisax=subplot(length(VOIs),length(r_cond),(vv-1)*length(condstr)+cc); hold on;
            plot(myth,thisrecon,'-','Color',mycolors(cc,:),'LineWidth',1);
            plot(myth,thisci,'--','Color',mycolors(cc,:),'LineWidth',0.5);
            
         
            
            if cc == 1
                
                ylabel(VOIs{vv});
            end
            
            if vv == 1
                title(condstr{cc});
            end
            
             set(thisax,'XTick',[-90 0 90],'TickLength',[0.025 0.025],'Box','off','XLim',[min(myth) max(myth)]);
            
            if cc~=1
                set(thisax,'YTickLabel',[]);
            end
            
            if vv~=length(VOIs)
                set(thisax,'XTickLabel',[]);
            end
            
 
            axi(end+1) = thisax; clear thisax thisidx thisrecon thisbf thisbf_params;
        end
        
        
        set(gcf,'Position',[440   365   463   433],'NumberTitle','off','Name',sprintf('%s',delay_str{tt}));
        
 
        
    end
end


set(axi,'YLim',[-0.1 0.3]);



return