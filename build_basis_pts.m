% build_basis_pts.m

% builds basis set given just a list of rf points, basis fcn sizes, and x,y grid along which to evaluate
% (doesn't assume grid of RFs, etc)
%
% built for wmDrop - hexagonal basis set
%
% allows for mutliple filter sizes, could use for multiscale RF estimation
%
% assumes we're using make2dcos....
%
% just a wrapper around /gridfit/make_grid.m
%
% Tommy Sprague (TCS), 10/24/14



function basis_set = build_basis_pts(centersX,centersY,basis_size,xx,yy)



%root = load_root;
%addpath([root 'mFiles/gridfit']);

if length(basis_size)==1
    basis_size = basis_size*ones(length(centersX),1);
end

basis_set = make_grid(@make2dcos_grid,[xx yy],[centersX centersY basis_size]);


return