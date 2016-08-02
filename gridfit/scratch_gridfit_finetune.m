% scratch_gridfit_finetune.m

fn_fits = '/usr/local/serenceslab/tommy/encodeContrast/encodeContrast_pRFfits/AA61_PRFRidge_CrossCondMeanBIC_moreBasisFcns_GridFit1_V1';

load(fn_fits);
load(prffile);

xx = linspace(-fov(1)/2,fov(1)/2,res(1));
yy = linspace(-fov(2)/2,fov(2)/2,res(2));

[mx, my] = meshgrid(xx,yy);clear xx yy;

evalpts = [reshape(mx,numel(mx),1) reshape(my,numel(my),1)];

%parpool(8);

[bf_ft,ft_err,ft_fcn] = gridfit_finetune(prf_vec{1}.',@make2dcos_grid,bf{1},evalpts);

% convert size into FWHM
bf_ft(:,3) = 0.5.*rad2fwhm(bf_ft(:,3));

plot_all_prfs(prf_vec{1}.',bf_ft,evalpts);