
% EXAMPLE OF INCOMPRESSIBLE MATERIAL
% K Bulk modulus = 30 * A10
% Neohookean Material

fprintf('Adding to the path GeoPDES folders\n')
addpath(genpath('C:/Program Files/MATLAB/R2017b/geopdes'))
addpath(genpath('C:/Program Files/MATLAB/R2017b/nurbs'))
addpath(genpath('C:/Users/utente/desktop/TESI/NonLinear Elasticity_v2'))
fprintf('Done \n')


fprintf('######################## \n')
fprintf('Test for incompressibility\n')
fprintf('######################## \n')

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
problem_data.L =1;
L=10;
H = 1;
problem_data.geo_name = nrb4surf([0 0], [L 0], [0 H],  [L H]);

% Type of boundary conditions
% problem_data.nmnn_sides   = [];
% problem_data.drchlt_sides = [];

% Physical parameters
A10 = 80;
A01 = 0; %NEOHOOKEAN MATERIAL        
K = 2400; %Bulk modulus
problem_data.mat_property = [A10, A01, K];
problem_data.lambda_x = 99.9/100;
 
% Source and boundary terms
    %Dirichlet B.C.
      problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
%      hx = @(x, y, ind) (0*x);
%      hy = @(x, y, ind) (0*x);
%      problem_data.h = @(x, y, ind) cat(1, ...
%                  reshape (hx (x,y, ind), [1, size(x)]), ...
%                  reshape (hy (x,y, ind), [1, size(x)]));
     
     %LATO DX
     problem_data.u1 = (1-problem_data.lambda_x)*L/2;
     %LATO SX
     problem_data.u2 = -(1-problem_data.lambda_x)*L/2;
    
%     %Neumann b.c.
%     gx = @(x, y, ind) (0*x);
%     gy = @(x, y, ind) (0*x-10);
%     problem_data.g = @(x, y, ind) cat(1, ...
%                 reshape (gx (x,y, ind), [1, size(x)]), ...
%                 reshape (gy (x,y, ind), [1, size(x)]));
% 
%      gx = @(x, y, ind) (0*x+lambda);
%      gy = @(x, y, ind) (0*x);
%      problem_data.gsx = @(x, y, ind) cat(1, ...
%                  reshape (gx (x,y, ind), [1, size(x)]), ...
%                  reshape (gy (x,y, ind), [1, size(x)]));
%      gx = @(x, y, ind) (0*x-lambda);
%      gy = @(x, y, ind) (0*x);
%      problem_data.gdx = @(x, y, ind) cat(1, ...
%                  reshape (gx (x,y, ind), [1, size(x)]), ...
%                  reshape (gy (x,y, ind), [1, size(x)]));
%     
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
method_data.nsub       = [60 60];     % Number of subdivisions
method_data.nquad      = [2 2];     % Points for the Gaussian quadrature rule
method_data.nel_small =0;


method_data.eps_d = 1e-4;
method_data.eps_r = 1e-4;
method_data.num_max_it = 100;


% 3) CALL TO THE SOLVER
fprintf('######################## \n')
fprintf('Solving NON linear system \n')
fprintf('######################## \n')
[geometry, msh, space, u, errore] = solve_NON_linear_elasticity (problem_data, method_data);
fprintf('######################## \n')
fprintf('NON linear system solved\n')
fprintf('######################## \n \n')

% 4) POST-PROCESSING
vtk_pts = {linspace(0, 1, 100), linspace(0, 1, 100)}; 

output_file = 'Displacement';
output_file_grad = 'Displ_grad';
fprintf('###################################### \n')
fprintf ('Displacement being saved in: %s.vts\n', output_file)

%Vtk of displacement
sp_to_vtk (u, space, geometry, vtk_pts, output_file, {'displacement'}, {'value'})

%Vtk of gradient
fprintf ('Gradient being saved in: %s.vts\n', output_file_grad)
sp_to_vtk (u, space, geometry, vtk_pts, output_file_grad, {'gradient'}, {'gradient'})

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
[eu, F] = sp_eval (u, space, geometry, vtk_pts, {'value', 'gradient'});
[sigma_stress, S] = stress_eval (u, space, geometry, vtk_pts, problem_data.mat_property);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

second_piola_file_1= 'Sxx.vtk';
second_piola_file_2= 'Sxy.vtk';
second_piola_file_3= 'Syy.vtk';
fprintf ('2-nd Piola Stress tensor being saved in: %s, %s, %s \n', second_piola_file_1,second_piola_file_2,second_piola_file_3)
vtkwrite(second_piola_file_1, 'structured_points', 'S_{xx}', squeeze(S(1,1,:,:)))
vtkwrite(second_piola_file_2, 'structured_points', 'S_{xy}', squeeze(S(1,2,:,:)))
vtkwrite(second_piola_file_3, 'structured_points', 'S_{yy}', squeeze(S(2,2,:,:)))

cauchy_file_1= 'Sigmaxx.vtk';
cauchy_file_2= 'Sigmaxy.vtk';
cauchy_file_3= 'Sigmayy.vtk';
fprintf ('Cauchy Stress tensor being saved in: %s, %s, %s \n', cauchy_file_1,cauchy_file_2,cauchy_file_3)
vtkwrite(cauchy_file_1, 'structured_points', 'sigma_{xx}', squeeze(sigma_stress(1,1,:,:)))
vtkwrite(cauchy_file_2, 'structured_points', 'sigma_{xy}', squeeze(sigma_stress(1,2,:,:)))
vtkwrite(cauchy_file_3, 'structured_points', 'sigma_{yy}', squeeze(sigma_stress(2,2,:,:)))
fprintf('###################################### \n')
 
figure
subplot(2,2,1)
surf (X, Y, squeeze(eu{2}(1,1,:,:)))
view(2); shading interp; colorbar;
title ('Gradient component \nabla_{xx}')

subplot(2,2,2)
surf (X, Y, squeeze(eu{2}(1,2,:,:)))
view(2); shading interp; colorbar;
title ('Gradient component \nabla_{xy}')

subplot(2,2,3)
surf (X, Y, squeeze(eu{2}(2,1,:,:)))
view(2); shading interp; colorbar;
title ('Gradient component \nabla_{yx}')

subplot(2,2,4)
surf (X, Y, squeeze(eu{2}(2,2,:,:)))
view(2); shading interp; colorbar;
title ('Gradient component \nabla_{yy}')

figure
subplot(1,3,1)
surf (X, Y, squeeze(S(1,1,:,:)))
view(2); shading interp; colorbar;
title ('2-nd Piola Kirchhoff stress component S_{xx}')

subplot(1,3,2)
surf (X, Y, squeeze(S(1,2,:,:)))
view(2); shading interp; colorbar;
title ('2-nd Piola Kirchhoff stress component S_{xy}')


subplot(1,3,3)
surf (X, Y, squeeze(S(2,2,:,:)))
view(2); shading interp; colorbar;
title ('2-nd Piola Kirchhoff stress component S_{yy}')

figure
subplot(1,3,1)
surf (X, Y, squeeze(sigma_stress(1,1,:,:)))
view(2); shading interp; colorbar;
title ('Cauchy stress component \sigma_{xx}')

subplot(1,3,2)
surf (X, Y, squeeze(sigma_stress(1,2,:,:)))
view(2); shading interp; colorbar;
title ('Cauchy stress component \sigma_{xy}')


subplot(1,3,3)
surf (X, Y, squeeze(sigma_stress(2,2,:,:)))
view(2); shading interp; colorbar;
title ('Cauchy stress component \sigma_{yy}')
