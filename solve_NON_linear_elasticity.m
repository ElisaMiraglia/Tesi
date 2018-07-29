
% SOLVE_NON_LINEAR_ELASTICITY: Solve a NON linear elasticity problem on a NURBS domain.
%
% The function solves the linear elasticity problem
%
%      - div (F*S(u)) = f       in Omega = (0,1)^2
%      F*S(u) \cdot n = g       on Gamma_N
%                   u = h       on Gamma_D
%
%   u:            displacement vector
%   F:            deformation gradient
%   S:            second Piola-Kirchhoff stress tensor
%
% USAGE:
%
%  [geometry, msh, space, u] = solve_linear_NON_elasticity (problem_data, method_data)
%
% INPUT:
%
%  problem_data: a structure with data of the problem. It contains the fields:
%    - geo_name:     name of the file containing the geometry
%    - nmnn_sides:   sides with Neumann boundary condition (may be empty)
%    - drchlt_sides: sides with Dirichlet boundary condition
%    - f:            source term
%    - h:            function for Dirichlet boundary condition
%    - g:            function for Neumann condition (if nmnn_sides is not empty)
%
%  method_data : a structure with discretization data. Its fields are:
%    - degree:     degree of the spline functions.
%    - regularity: continuity of the spline functions.
%    - nsub:       number of subelements with respect to the geometry mesh 
%                   (nsub=1 leaves the mesh unchanged)
%    - nquad:      number of points for Gaussian quadrature rule
%    - eps_d:      tollerance on solution of Newton iteration
%    - eps_r:      tollerance on residual
%    - num_max_it  Max iterations of Newton Algorithm
%
% OUTPUT:
%
%  geometry: geometry structure (see geo_load)
%  msh:      mesh object that defines the quadrature rule (see msh_cartesian)
%  space:    space object that defines the discrete basis functions (see sp_vector)
%  u:        the computed degrees of freedom
%

function [geometry, msh, sp, u] = ...
              solve_NON_linear_elasticity (problem_data, method_data)

% Extract the fields from the data structures into local variables
data_names = fieldnames (problem_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= problem_data.(data_names{iopt});']);
end
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
  eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end


% % Construct geometry structure
geometry = geo_load (geo_name);
degelev  = max (degree - (geometry.nurbs.order-1), 0);
nurbs    = nrbdegelev (geometry.nurbs, degelev);
[~, ~, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

nurbs    = nrbkntins (nurbs, nknots);
geometry = geo_load (nurbs);


% Construct msh structure
rule     = msh_gauss_nodes (nquad);
[qn, qw] = msh_set_quad_nodes (geometry.nurbs.knots, rule);
msh      = msh_cartesian (geometry.nurbs.knots, qn, qw, geometry);


% Construct space structure
space_scalar = sp_nurbs (nurbs, msh);
scalar_spaces = repmat ({space_scalar}, 1, msh.rdim);
sp = sp_vector (scalar_spaces, msh);
clear space_scalar scalar_spaces



%Vector for evaluation of the matrix
x_1 = [];
x_2 = [];
     for i = 1:size(msh.qn{1},2)  
        x_1 = [x_1, msh.qn{1}(:,i)'];
     end

    for j = 1:size(msh.qn{2},2)  
        x_2 = [x_2, msh.qn{2}(:,j)'];
    end
x = {x_1, x_2};

%Initial guess
u = zeros (sp.ndof, 1);

num_col = method_data.nquad(1);     % Points for the Gaussian quadrature rule along x direction
num_row = method_data.nquad(2);     % Points for the Gaussian quadrature rule along y direction


    % Assemble the rhs which is constant (ext force + Neumann condition)
    b    = op_f_v_tp (sp, msh, f);
   
    % Apply Neumann boundary conditions   
%   iside=1;
%   g_sx_side = @(varargin) gsx(varargin{:},iside);
%   dofs = sp.boundary(iside).dofs;
%   b(dofs) = b(dofs) + op_f_v_tp (sp.boundary(iside), msh.boundary(iside), g_sx_side);
% 
%   iside=2;
%   g_dx_side = @(varargin) gdx(varargin{:},iside);
%   dofs = sp.boundary(iside).dofs;
%   b(dofs) = b(dofs) + op_f_v_tp (sp.boundary(iside), msh.boundary(iside), g_dx_side);

  lambda_x = 95/100;
    
    %Apply Dirichlet boundary conditions
    [u_drchlt, drchlt_dofs1] = sp_drchlt_l2_proj (sp, msh, h, [1]);
    %drchlt_dofs1
    dim = size(drchlt_dofs1,1);
    u(drchlt_dofs1(1:dim/2)) = (1-lambda_x)*L/2;
    
    drchlt_dofs = [drchlt_dofs1(1:dim/2)]';
    [u_drchlt, drchlt_dofs2] = sp_drchlt_l2_proj (sp, msh, h, [2]);
    dim = size(drchlt_dofs2,1);
    u(drchlt_dofs2(1:dim/2)) = -(1-lambda_x)*L/2;
    drchlt_dofs = [drchlt_dofs, drchlt_dofs2(1:dim/2)'];
    
    [u_drchlt, drchlt_dofs3] = sp_drchlt_l2_proj (sp, msh, h, [3]);
    dim = size(drchlt_dofs3,1);
    u(drchlt_dofs3(dim/2+1:dim)) = 0;
    
    drchlt_dofs = [drchlt_dofs, drchlt_dofs3(dim/2+1:dim)'];
    sort(drchlt_dofs);
    u_drchlt = zeros(size(drchlt_dofs,1));
    %[u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp, msh, h, drchlt_sides);
    
    %u(drchlt_dofs) = u_drchlt;
       
    %Initialization of parameters for Newton cycle    
    err_d = 1000;
    err_r = 1000;
    num_it=0;
    rhs = zeros(size(b));
    
    %mat_property = [mat_property2,mat_property1]; %small film, big film
    
    %N-R cycle
    while err_d > eps_d && err_r > eps_r && num_it<num_max_it
             [val, ~] = sp_eval (u, sp, geometry, x, {'value', 'gradient'});
              D = val{2};

             num_it=num_it+1;
             fprintf('It N-R: %d \n', num_it);
             
             % Apply Dirichlet boundary conditions, for variation they are
             % always 0
             delta_u = zeros (sp.ndof, 1);            
             delta_u(drchlt_dofs) = 0*u_drchlt;
             int_dofs = setdiff (1:sp.ndof, [drchlt_dofs]);
             
             %   int_dofs = setdiff (1:sp.ndof, [drchlt_dofs((dim/2+1):dim)]);

             mat    = op_mat_stiff_tp (sp, sp, msh, D, num_row, num_col, mat_property, nel_small)+op_geo_stiff_tp (sp, sp, msh, D, num_row, num_col, mat_property, nel_small);

             f_s  = op_f_d_s_tp(sp, msh, D, num_row, num_col, mat_property, nel_small);             
             rhs(int_dofs)  = f_s(int_dofs) -b(int_dofs); 
             
             % Solve the linearyzed system
             delta_u(int_dofs) = - mat(int_dofs, int_dofs) \ rhs(int_dofs);
             err_d = norm(delta_u)/norm(u)
             err_r = norm(rhs)/norm(b)
        
             %update the solution
             u = u + delta_u;
             
       
%               vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)}; 
%               [eu, F] = sp_eval (u, sp, geometry, vtk_pts, {'value'});
%               [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
% 
%               quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
%               axis equal tight
%               title (['Computed displacement IT: ', num2str(num_it)])
%               %pause
    
    end
    %fprintf("Newton Cycle stopped with err_d %f and err_r %f \n", err_d, err_r);
end