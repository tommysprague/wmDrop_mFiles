% gridfit_noregress.m
% TCS 2/9/2015, adapted from gridfit.m
function [bf_params, err, bf_fcn] = gridfit_noregress(data,grid,grid_params,use_gpu)

% data should be datapts x things to fit (e.g., voxels)
% grid should be datapts x predictors
% grid_params is the params used to make grid - just used to output
% bf_params as the last two params, and used in case
% optional fminsearch optimization step at end is used (?)
%
% uses parfor, but assumes the parallel cluster has been initialized
% already in parent funciton
%
% does not compute a/b regression - just computes R2/error for all
% predictors


if nargin < 4
    use_gpu = 0;
end

%myones = ones(size(grid,1),1);
%allcoeffs = nan(size(grid,2),2,size(data,2));
allerr = nan(size(grid,2),size(data,2));


%convert things to gpuarrays
if use_gpu==1
    data = gpuArray(data);
    grid = gpuArray(grid);
    %myones = gpuArray(myones);
    %allcoeffs = gpuArray(allcoeffs);%nan(size(grid,2),size(data,2),2));
    allerr = gpuArray(allerr);%nan(size(grid,2),size(data,2))); % grid x vox
end


% data should be datapts x things to fit (e.g., voxels)
% grid should be datapts x predictors

% alldata will be datapts x thingstofit x n_predictors (repeated)
alldata = repmat(data,1,1,size(predictors,2));
allgrid = repmat(permute(grid,[size(grid,1),1,size(grid,2)]),1,size(data,2),1);



% prf_vec{1} is n_vox x n_datapts
% grid is n_datapts x n_predictors
tic;
parfor ii = 1:size(grid,2)
%for ii = 1:size(grid,2)
    fprintf('grid iter %i\n',ii);
    
    %allcoeffs(ii,:,:) = [grid(:,ii) myones]\data; % returns predictors(2) x nvox 
    
    %allerr(ii,:) = sqrt(mean((data - [grid(:,ii) myones]*squeeze(allcoeffs(ii,:,:))).^2,1));
    allerr(ii,:) = sqrt(mean((data - [repmat(grid(:,ii),1,size(data,2))]).^2,1));
    
end
toc;
%allamp = squeeze(allcoeffs(:,1,:));
%allerr(allamp < amp_range(1) | allamp > amp_range(2)) = inf;
% now find the best fit and arrange those values for returning
[err,fidx] = min(allerr,[],1); % 1 x vox

if use_gpu==1
%    allcoeffs = gather(allcoeffs);
    grid = gather(grid);
    err = gather(err);
end

bf_params = nan(size(data,2),size(grid_params,2)+2);
bf_fcn = nan(size(data,2),size(grid,1)); % vox x datapts
for ii = 1:length(fidx) % for each voxel
    bf_params(ii,:) = [grid_params(fidx(ii),:)];
    bf_fcn(ii,:)  = grid(:,fidx(ii))';
end




%% TODO: fine-tune fits w/ fminsearch


return