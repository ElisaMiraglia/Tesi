%Plane strain compression with mu_f/mu_s=30 and no substrate prestretch
% Neohookean Material

fprintf('Adding to the path GeoPDES folders\n')
addpath(genpath('C:/Program Files/MATLAB/R2017b/geopdes'))
addpath(genpath('C:/Program Files/MATLAB/R2017b/nurbs'))
addpath(genpath('C:/Users/utente/desktop/TESI/NonLinear Elasticity'))
fprintf('Done \n')


fprintf('######################## \n')
fprintf('Test with 2 films\n')
fprintf('######################## \n')

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
problem_data.small_film = 2;
problem_data.big_film = 62;
H = problem_data.small_film + problem_data.big_film;
L = 4*H;
problem_data.geo_name = nrb4surf([0 0], [L 0], [0 H],  [L H]);


% Type of boundary conditions
problem_data.nmnn_sides   = [];
problem_data.drchlt_sides = [3];

% Physical parameters
%Big film
A10 = 80;
A01 = 0; %NEOHOOKEAN MATERIAL        
K = 2400; %Bulk modulus (30*A10)
problem_data.mat_property2 = [A10, A01, K];
%Small film
problem_data.mat_property1 = 30.*problem_data.mat_property2;
problem_data.mat_property = [problem_data.mat_property1,problem_data.mat_property2];
% Source and boundary terms
    %Dirichlet B.C.
    %problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
     hx = @(x, y, ind) (0*x);
     hy = @(x, y, ind) (0*x);
     problem_data.h = @(x, y, ind) cat(1, ...
                 reshape (hx (x,y, ind), [1, size(x)]), ...
                 reshape (hy (x,y, ind), [1, size(x)]));

    %Neumann b.c.
    gx = @(x, y, ind) (0*x);
    gy = @(x, y, ind) (0*x-10);
    problem_data.g = @(x, y, ind) cat(1, ...
                reshape (gx (x,y, ind), [1, size(x)]), ...
                reshape (gy (x,y, ind), [1, size(x)]));

     lambda = 50;
     gx = @(x, y, ind) (0*x+lambda);
     gy = @(x, y, ind) (0*x);
     problem_data.gsx = @(x, y, ind) cat(1, ...
                 reshape (gx (x,y, ind), [1, size(x)]), ...
                 reshape (gy (x,y, ind), [1, size(x)]));
     gx = @(x, y, ind) (0*x-lambda);
     gy = @(x, y, ind) (0*x);
     problem_data.gdx = @(x, y, ind) cat(1, ...
                 reshape (gx (x,y, ind), [1, size(x)]), ...
                 reshape (gy (x,y, ind), [1, size(x)]));
    
     %Source term
     fx = @(x, y) (0*x);
     fy = fx;
     problem_data.f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));
            
            
% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [1 1];     % Degree of the bsplines
method_data.regularity = [0 0];     % Regularity of the splines
method_data.nsub       = [62 62];     % Number of subdivisions
method_data.nquad      = [2 2];     % Points for the Gaussian quadrature rule

method_data.nel_small = floor(H/method_data.nsub(1)*problem_data.small_film); 

method_data.eps_d = 1e-4;
method_data.eps_r = 1e-4;
method_data.num_max_it = 20;


% 3) CALL TO THE SOLVER
fprintf('######################## \n')
fprintf('Solving NON linear system \n')
fprintf('######################## \n')
[geometry, msh, space, u] = solve_NON_linear_elasticity (problem_data, method_data);
fprintf('######################## \n')
fprintf('NON linear system solved\n')
fprintf('######################## \n \n')

% 4) POST-PROCESSING
output_file = 'Plane_strain_compression_mu_rate_30';
fprintf('###################################### \n')
fprintf ('Results being saved in: %s \n', output_file)
fprintf('###################################### \n')

%vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
vtk_pts = {linspace(0, L, method_data.nsub(1)), linspace(0, H, method_data.nsub(2))};
sp_to_vtk (u, space, geometry, vtk_pts, output_file, {'displacement'}, {'value'}, {'gradient'})

figure
subplot (1, 2, 1)
sp_plot_solution (u, space, geometry, vtk_pts)
axis equal tight
title ('Computed solution')
subplot (1, 2, 2)
def_geom = geo_deform (u, space, geometry);
nrbplot (def_geom.nurbs, method_data.nsub)
view(2)
title ('Deformed configuration')

%stress computation
stress_pts = {linspace(0, L, method_data.nsub(1)), linspace(0, H, method_data.nsub(2))};
[eu, F] = sp_eval (u, space, geometry, stress_pts, {'value', 'gradient'});
S = stress_eval_2film (u, space, geometry, stress_pts, problem_data.mat_property);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));


figure
surf (X, Y, squeeze(eu{2}(1,1,:,:)))
view(2); axis equal; shading interp; colorbar;
title ('Gradient component \nabla_{xx}')

figure
surf (X, Y, squeeze(S(1,1,:,:)))
view(2); axis equal; shading interp; colorbar;
title ('Nominal stress component S_{xx}')
