% EXAMPLE

% 1) PHYSICAL DATA OF THE PROBLEm
clear problem_data
L = 1;
H = 1;
problem_data.geo_name = nrb4surf([0 0], [L 0], [0 H],  [L H]);

% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [1 2 3 4];
problem_data.press_sides = [];
problem_data.symm_sides = [];

% Physical parameters (Mooney material)
A10 = 8000; %Ordine di 10^3 (controllare)
A01 = 2000;
        
K = 0; %Bulk modulus

% Physical parameters
E  =  1; nu = .3; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));


problem_data.mat_property = [A10, A01, K];

% Source and boundary terms
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2)); %Dirichlet b.c.
problem_data.g = @(x, y, ind) ones (2, size (x, 1), size (x, 2)); %Neumann b.c.

%fx = @(x, y) -(-(problem_data.mu_lame(x, y)*3 + problem_data.lambda_lame(x, y)).*sin(2*pi*x).*sin(2*pi*y) + ...
 %    (problem_data.mu_lame(x, y) + problem_data.lambda_lame(x, y)).*cos(2*pi*x).*cos(2*pi*y))*(2*pi)^2;
 fx = @(x, y) (0*x-1);
 fy = @(x, y) (0*x -1);
problem_data.f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));
            
%problem_data.h       = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));



% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
%method_data.degree     = [3 3];     % Degree of the bsplines
%method_data.regularity = [2 2];     % Regularity of the splines
%method_data.nsub       = [9 9];     % Number of subdivisions
%method_data.nquad      = [4 4];     % Points for the Gaussian quadrature rule
method_data.degree     = [1 1];     % Degree of the bsplines
method_data.regularity = [0 0];     % Regularity of the splines
method_data.nsub       = [4 4];     % Number of subdivisions
method_data.nquad      = [2 2];     % Points for the Gaussian quadrature rule

method_data.eps_d = 0.001;
method_data.eps_r = 0.001;
method_data.num_max_it = 5;


% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_NON_linear_elasticity (problem_data, method_data);
%[geometry, msh, space, u] = solve_linear_elasticity (problem_data, method_data);

% % 4) POST-PROCESSING.
 vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)}; 
 [eu, F] = sp_eval (u, space, geometry, vtk_pts, {'value'});
 [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

figure
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal tight
title ('Computed displacement')



%def_geom = geo_deform (u, space, geometry);
%figure
%nrbplot (def_geom.nurbs, method_data.nsub, 'light', 'on')
%axis equal tight
%title ('Deformed configuration')
