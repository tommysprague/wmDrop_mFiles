% gridfit_scratch.m

% playing around with using make_grid.m and make2dcos_grid.m




xx = linspace(-6,6,121);
yy = linspace(-3,3,61);

[mx, my] = meshgrid(xx,yy);

evalpts = [reshape(mx,numel(mx),1) reshape(my,numel(my),1)];

% xgrid = linspace(-6,6,25);
% ygrid = linspace(-3,3,13);
% sgrid = linspace(1.8,8,20);

xgrid = linspace(-6,6,13);
ygrid = linspace(-3,3,7);
sgrid = linspace(1.8,8,10);


grid_params =  nan(length(xgrid)*length(ygrid)*length(sgrid),3);

idx = 1;
for xx = 1:length(xgrid)
    for yy = 1:length(ygrid)
        for ss = 1:length(sgrid)
            grid_params(idx,:) = [xgrid(xx) ygrid(yy) sgrid(ss)];
            idx = idx+1;
        end
    end
end

mygrid = make_grid(@make2dcos_grid,evalpts,grid_params);

load('/usr/local/serenceslab/tommy/encodeContrast/encodeContrast_pRFs/AA61_PRFRidge_CrossCondMeanBIC_moreBasisFcns_V1.mat');
data = prf_vec{2}';
%mygrid = 
[a,b,c] = gridfit(data,mygrid,grid_params,[0 inf]);