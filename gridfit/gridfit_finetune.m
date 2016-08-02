% gridfit_finetune.m

function [bf_ft,bf_err,bf_fcn,ex_flag] = gridfit_finetune(data,fitfcn,bf_grid,evalpts)

% refines fits computed with gridfit.m using fminsearch and fitfcn, which
% is used here to generate an error function (root mean squared error)
%
% data formatted like gridfit.m:
% - data: n_datapts x n_things_to_fit (voxels)
% - fitfcn: function handle to a normalized-amplitude/baseline function
%   prototype used for GLM fitting in gridfit fitting (base = 0, amp = 1)
% - bf_grid: the best-fit from gridfit.m (its return argument bf_params)
% - evalpts: same as gridfit.m, the points at which fitfcn is evaluated for
%   a given set of parameters (e.g, x/y coords) - n_datapts x n_coords

% assumes bf_grid is output from gridfit.m - last 2 columns are amplitude
% and baseline, respectively
% TODO: optional argument that specifies the column for amp, column for
% baseline, and defaults to just end-1, end (current implementation does
% this)



% TODO: some basic checks to be sure the data is all formatted correctly





% initialize all variables before parallel
bf_ft = nan(size(bf_grid));
bf_err = nan(size(bf_grid,1),1);
bf_fcn = nan(size(data));
ex_flag = nan(size(bf_grid,1),1);

amp_idx = size(bf_grid,2)-1;
base_idx = size(bf_grid,2);

parfor vv = 1:size(data,2)

    
    if mod(vv,10)==0
        fprintf('Voxel %i\n',vv);
    end
    
    d = data(:,vv);
    
    % here's the error function (RMSE)
    err_fcn = @(p) sqrt(mean( (   (p(end-1)*fitfcn(evalpts,p)+p(end)) - d).^2 ) );

    
    [bf_ft(vv,:),bf_err(vv),ex_flag(vv)] = fminsearch(err_fcn,bf_grid(vv,:));
    this_fits = bf_ft(vv,:);
    this_amp = this_fits(amp_idx);
    this_base = this_fits(base_idx);
    bf_fcn(:,vv) = this_amp*fitfcn(evalpts,bf_ft(vv,:)) + this_base;
    
end



return