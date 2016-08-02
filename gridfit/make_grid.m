function grid = make_grid(gridfcn,evalpts,grid_params)

% gridfcn is an anonymous function used to generate the surfaces - it takes
% two inputs, evalpts (n_pts x n_coords, like x, y coords over a surface)
% and grid_params(ii,:) and returns a n_pts x 1 column vector corresponding
% to gridfcn evaluated at evalpts with grid_params(ii,:) inputs
%
% evalpts is the set of points (each row is a point) at which to evaluate
% gridfcn for each set of parameter values (e.g., x and y coords)
%
% grid_params is a n_total_fits x n_params matrix, where each row is the set of
% inputs to gridfcn for that particular surface (evaluated at evalpts)


grid = nan(size(evalpts,1),size(grid_params,1));


for ii = 1:size(grid_params,1)
    %dd = gridfcn(evalpts,grid_params(ii,:));
    grid(:,ii) = gridfcn(evalpts,grid_params(ii,:));
end


return