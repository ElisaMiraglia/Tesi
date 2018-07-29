
% EXAMPLE

fprintf('Adding to the path GeoPDES folders\n')
addpath(genpath('C:/Program Files/MATLAB/R2017b/geopdes'))
addpath(genpath('C:/Program Files/MATLAB/R2017b/nurbs'))
addpath(genpath('C:/Users/utente/desktop/TESI/NonLinear Elasticity'))
fprintf('Done \n')


fprintf('######################## \n')
fprintf('Back test with Linear problem\n')
fprintf('######################## \n')

% 1) PHYSICAL DATA OF THE PROBLEM
clear problem_data
L = 1;
H = 1;
problem_data.geo_name = nrb4surf([0 0], [L 0], [0 H],  [L H]);

% Type of boundary conditions
problem_data.nmnn_sides   = [4];
problem_data.drchlt_sides = [1];
problem_data.press_sides = [];
problem_data.symm_sides = [];

% Physical parameters (Mooney material)
A10 = 80; %Ordine di 10^3 (controllare)
A01 = 20;
        
K = 0; %Bulk modulus

% Physical parameters
E  =  1; nu = .3; 
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));


problem_data.mat_property = [A10, A01, K];

% Source and boundary terms
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2)); %Dirichlet b.c.
%problem_data.g = @(x, y, ind) zeros (2, size (x, 1), size (x, 2)); %Neumann b.c.
gx = @(x, y, ind) (0*x);
gy = @(x, y, ind) (0*x-10);
problem_data.g = @(x, y, ind) cat(1, ...
                reshape (gx (x,y, ind), [1, size(x)]), ...
                reshape (gy (x,y, ind), [1, size(x)]));

% fx = @(x, y) -(-(problem_data.mu_lame(x, y)*3 + problem_data.lambda_lame(x, y)).*sin(2*pi*x).*sin(2*pi*y) + ...
%      (problem_data.mu_lame(x, y) + problem_data.lambda_lame(x, y)).*cos(2*pi*x).*cos(2*pi*y))*(2*pi)^2;
% fy = fx;
%fx = @(x, y) (sin(2*pi*x).*sin(2*pi*y));
%fy = @(x, y) (cos(2*pi*x).*cos(2*pi*y));


fx = @(x, y) (0*x);
fy = fx;
problem_data.f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));
            
%problem_data.h       = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));



% 2) CHOICE OF THE DISCRETIZATION PARAMETERS
clear method_data
method_data.degree     = [1 1];     % Degree of the bsplines
method_data.regularity = [0 0];     % Regularity of the splines
method_data.nsub       = [10 10];     % Number of subdivisions
method_data.nquad      = [2 2];     % Points for the Gaussian quadrature rule

method_data.eps_d = 1e-5;
method_data.eps_r = 1e-5;
method_data.num_max_it = 100;


% 3) CALL TO THE SOLVER

fprintf('######################## \n')
fprintf('Solving first iteration of NON linear system \n')
fprintf('######################## \n')
[geometry, msh, space, u] = solve_NON_linear_elasticity (problem_data, method_data);
fprintf('NON linear system solved\n')

fprintf('######################## \n')
fprintf('Solving linear system \n')
fprintf('######################## \n')
[geometrylin, mshlin, spacelin, ulin] = solve_linear_elasticity (problem_data, method_data);

% % 4) POST-PROCESSING.
 vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)};
 [eu, F] = sp_eval (u, space, geometry, vtk_pts, {'value'});
 [eulin, Flin] = sp_eval (ulin, space, geometry, vtk_pts, {'value'});
 [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
 [Xlin, Ylin]  = deal (squeeze(Flin(1,:,:)), squeeze(Flin(2,:,:)));

 
figure
subplot(1,3,1)
quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
axis equal tight
title ('Computed displacement NON linear solver')
subplot(1,3,2)
quiver (Xlin, Ylin, squeeze(eulin(1,:,:)), squeeze(eulin(2,:,:)))
axis equal tight
title ('Computed displacement linear solver')
err = sqrt((eulin(1,:,:)-eu(1,:,:)).^2 +(eulin(2,:,:)-eu(2,:,:)).^2);
err = reshape(err,size(vtk_pts{1},2),size(vtk_pts{1},2));
subplot(1,3,3)
mesh(X, Y, err)
colorbar('southoutside')
view(90,90)
axis equal tight
title ('Error')

fprintf('Error: u-ulin %f \n', norm(u-ulin))