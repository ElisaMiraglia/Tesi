%BOUSSINESQ ON RECTANGULAR DOMAIN [0,L]x[0,H]

%1) PHYSICAL DATA
clear problem_data
nx = 10;
ny = 6;
L = 3;
H = 3;
eps = L./nx;
problem_data.geo_name = nrb4surf([-L H], [L H], [-L 0],  [L 0]);

problem_data.nmnn_sides = []; %Neumann
problem_data.press_sides  = [3]; %Pressure
problem_data.drchlt_sides = [1 2 4]; %Dir
problem_data.symm_sides   = []; %Symmetry

E=1;
nu = .3;
problem_data.lambda_lame = @(x, y) ((nu*E)/((1+nu)*(1-2*nu)) * ones (size (x)));
problem_data.mu_lame = @(x, y) (E/(2*(1+nu)) * ones (size (x)));

problem_data.L = 3;
problem_data.H = 3;

fx = @(x, y) 0.*x;
fy = @(x, y) 0.*x;
problem_data.f = @(x, y) cat(1, ...
                reshape (fx (x,y), [1, size(x)]), ...
                reshape (fy (x,y), [1, size(x)]));

%Boundary terms:
problem_data.h = @(x, y, ind) zeros (2, size (x, 1), size (x, 2));
problem_data.g = @(x, y, ind) zeros(2, size (x, 1), size (x, 2));
%problem_data.p = @(x, y, ind) 1.* (x<=(0.1).*L).*(x>=(-0.1).*L) .* ones (size (x));
problem_data.p = @(x, y, ind) ones (size (x));


 % 2) CHOICE OF THE DISCRETIZATION PARAMETERS

method_data.degree       = [1  1];  % Degree of the splines (pressure space)
method_data.regularity   = [0  0];  % Regularity of the splines (pressure space)
method_data.nsub         = [1  1];  % Number of subdivisions
method_data.nquad        = [1  1];  % Points for the Gaussian quadrature rule

% 3) CALL TO THE SOLVER
[geometry, msh, space, u] = solve_linear_elasticity(problem_data, method_data);


% 4) POST-PROCESSING. 
% 4.1) Export to Paraview
output_file = 'plane_strain_square_Deg3_Reg2_Sub9';

vtk_pts = {linspace(0, 1, 21), linspace(0, 1, 21)};
fprintf ('results being saved in: %s \n \n', output_file)
sp_to_vtk (u, space, geometry, vtk_pts, output_file, {'displacement', 'stress'}, {'value', 'stress'}, ...
    problem_data.lambda_lame, problem_data.mu_lame)

% 4.2) Plot in Matlab. Comparison with the exact solution.
[eu, F] = sp_eval (u, space, geometry, vtk_pts);
[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
quiver(X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
title ('Numerical solution'), axis equal tight

def_geom = geo_deform (u, space, geometry);
subplot(1,2,1)
nrbplot (def_geom.nurbs, [nx ny], 'light', 'on')
k=1;
figure(k)
nrbkntplot(def_geom.nurbs)
title ('Deformed configuration')
view(2)
k=k+1;
subplot(1,2,2)
nrbkntplot(geometry.nurbs)
% 'colormap', 'white'
view(2)
title ('Mesh grid')



% %%
% def_geom = geo_deform (u, sp, geometry);
% %subplot(1,2,1)
% %nrbplot (def_geom.nurbs, [nx ny], 'light', 'on')
% k=2;
% figure(k)
% nrbkntplot(def_geom.nurbs)
% title ('Deformed configuration')
% view(2)
% % title ('Mesh grid')
% k=k+1;
% % subplot(1,2,2)
% % nrbkntplot(geometry.nurbs)
% % % 'colormap', 'white'
% % view(2)
% 
