function plot_all_prfs( data, bf, evalpts, which_vox )
%PLOT_ALL_PRFS plotting function - shows all reconstructed pRFs, overlaid
%with their best-fit circle
%   Takes in outputs from gridfit.m or gridfit_finetune.m:
%   - DATA: n_datapts x n_voxels
%   - BF:   best fits, n_voxels x >= 3 (first 3 dimensions are x, y, size)
%   - XGRID: points along which to plot DATA (n_vert_pix x n_horiz_pix)
%   - YGRID:
%
%   Assumes that DATA is arranged as is typical w/ MESHGRID then RESHAPE -
%   uses first and last elements of EVALPTS as [X,Y] of first and last
%   element of DATA for each PRF

if nargin < 4
    which_vox = 1:size(data,2);
    %which_vox = 1:40;
end


% TODO: dynamically choose these...

n_cols = 8;
n_rows = 7;

myclim = [min(data(:)) max(data(:))];
myxlim = [min(evalpts(:,1)) max(evalpts(:,1))];
myylim = [min(evalpts(:,2)) max(evalpts(:,2))];

resx = length(unique(evalpts(:,1)));
resy = length(unique(evalpts(:,2)));

% used for drawing pRF fit on reconstructed pRF image
myt = linspace(2*pi/1000,2*pi,1000);
myx = cos(myt);
myy = sin(myt);

sp_idx = 1;
for vv = 1:length(which_vox)
    
    dd = which_vox(vv);
    
    if sp_idx == 1
        figure;
    end
    
    subplot(n_rows,n_cols,sp_idx);
    hold on;
    
    this_prf = reshape(data(:,dd),resy,resx);
    
    imagesc(evalpts(:,1),evalpts(:,2),this_prf,myclim);
    axis equal tight;
    plot(bf(dd,3)*myx+bf(dd,1),bf(dd,3)*myy+bf(dd,2),'--','Color',[1 1 1]);
    plot(bf(dd,1),bf(dd,2),'+','Color',[1 1 1],'MarkerSize',2);
    
    
    set(gca,'Box','off','XLim',myxlim,'YLim',myylim,'XTick',[-3.6 0 3.6], 'YTick',[-1.5 0 1.5],'XTickLabel',[],'YTickLabel',[],'TickLength',[0.025 0.025]);
    hold off;
    
    sp_idx = sp_idx+1;
    if sp_idx > n_rows * n_cols
        sp_idx = 1;
    end
    
    clear this_prf dd;
    
end



end

