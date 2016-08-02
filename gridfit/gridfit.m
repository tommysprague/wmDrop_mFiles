% gridfit.m

function [bf_params, err, bf_fcn] = gridfit(data,grid,grid_params,amp_range,use_gpu,use_parfor)

% data should be datapts x things to fit (e.g., voxels)
% grid should be datapts x predictors
% grid_params is the params used to make grid - just used to output
% bf_params including [a b] as the last two params, and used in case
% optional fminsearch optimization step at end is used (?)
%
% uses parfor, but assumes the parallel cluster has been initialized
% already in parent funciton

if nargin < 4
    amp_range = [-inf inf];
end

if nargin < 5
    use_gpu = 0;
end

if nargin < 6
    use_parfor = 0;
end

myones = ones(size(grid,1),1);
allcoeffs = nan(size(grid,2),2,size(data,2));
allerr = nan(size(grid,2),size(data,2));


%convert things to gpuarrays
if use_gpu==1
    data = gpuArray(data);
    grid = gpuArray(grid);
    myones = gpuArray(myones);
    allcoeffs = gpuArray(allcoeffs);%nan(size(grid,2),size(data,2),2));
    allerr = gpuArray(allerr);%nan(size(grid,2),size(data,2))); % grid x vox
end


% prf_vec{1} is n_vox x n_datapts
% grid is n_datapts x n_predictors
tic;
if use_parfor==0
%parfor ii = 1:size(grid,2)
for ii = 1:size(grid,2)
    %fprintf('grid iter %i\n',ii);
    
    allcoeffs(ii,:,:) = [grid(:,ii) myones]\data; % returns predictors(2) x nvox 
    %allerr(ii,:) = sqrt(mean((data - [grid(:,ii) myones]*squeeze(allcoeffs(ii,:,:))).^2,1));
    allerr(ii,:) = sqrt(mean((data - [grid(:,ii) myones]*shiftdim(allcoeffs(ii,:,:),1)).^2,1));
    
end

else
    parfor ii = 1:size(grid,2)
        if mod(ii,1000)==0
        fprintf('grid iter %i\n',ii);
        end
        allcoeffs(ii,:,:) = [grid(:,ii) myones]\data; % returns predictors(2) x nvox
        %allerr(ii,:) = sqrt(mean((data - [grid(:,ii) myones]*squeeze(allcoeffs(ii,:,:))).^2,1));
        allerr(ii,:) = sqrt(mean((data - [grid(:,ii) myones]*shiftdim(allcoeffs(ii,:,:),1)).^2,1));
        
    end

    
end

toc;
allamp = squeeze(allcoeffs(:,1,:));
allerr(allamp < amp_range(1) | allamp > amp_range(2)) = inf;
% now find the best fit and arrange those values for returning
[err,fidx] = min(allerr,[],1); % 1 x vox

if use_gpu==1
    allcoeffs = gather(allcoeffs);
    grid = gather(grid);
    err = gather(err);
end

bf_params = nan(size(data,2),size(grid_params,2)+2);
bf_fcn = nan(size(data,2),size(grid,1)); % vox x datapts
for ii = 1:length(fidx) % for each voxel
    bf_params(ii,:) = [grid_params(fidx(ii),:) squeeze(allcoeffs(fidx(ii),:,ii))];
    bf_fcn(ii,:)  = bf_params(ii,end-1)*grid(:,fidx(ii))' + bf_params(ii,end);
end

% TCS - updated above 6/17/2015 - wasn't computing amp*fit+b



return