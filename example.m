% EXAMPLE

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
L = 3;
H = 3;
problem_data.geo_name = nrb4surf([-L H], [L H], [-L 0],  [L 0]);

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [];

% Physical parameters (Mooney material)
A10 = 1;
A01 = 0;
K = 0; %Bulk modulus

problem_data.mat_property = [A10, A01, K];

% Source and boundary terms
problem_data.f = @(x, y) zeros (2, size (x, 1), size (x, 2));
problem_data.g = @(x, y) zeros (2, size (x, 1), size (x, 2));
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));

% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [3 3];     % Degree of the basis functions
method_data.regularity = [2 2];     % Regularity of the basis functions
method_data.nsub       = [4 2];     % Number of subdivisions
method_data.nquad      = [2 2];     % Points for the Gaussian quadrature rule
method_data.eps_d = 0.001;
method_data.eps_r = 0.001;


% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_NON_linear_elasticity (problem_data, method_data);

% % 4) POST-PROCESSING.
% vtk_pts = {linspace(0, 1, method_data.nsub(1)), linspace(0, 1, method_data,nsub(2))};
% 
% [eu, F] = sp_eval (u, space, geometry, vtk_pts, {'value'});
% [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
% 
% figure
% quiver (X, Y, squeeze(eu{1}(1,:,:)), squeeze(eu{1}(2,:,:)))
% axis equal tight
% title ('Computed displacement')
