% SOLVE_NON_LINEAR_ELASTICITY: Solve a NON linear elasticity problem on a NURBS domain.
%
% The function solves the linear elasticity problem
%
%      - div (F*S(u)) = f       in Omega = F((0,1)^n)
%      F*S(u) \cdot n = g       on Gamma_N
%                     u = h     on Gamma_D
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
%    - press_sides:  sides with pressure boundary condition (may be empty)
%    - symm_sides:   sides with symmetry boundary condition (may be empty)
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
%    - eps_
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
[rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);

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




%Creo il vettore per valutare la matrice

x_1 = [];
x_2 = [];
     for i = 1:size(msh.qn{1},2)  
        x_1 = [x_1, msh.qn{1}(:,i)'];
     end

    for j = 1:size(msh.qn{2},2)  
        x_2 = [x_2, msh.qn{2}(:,j)'];
    end

x = {x_1, x_2};

u = zeros (sp.ndof, 1);

num_col = method_data.nquad(1);     % Points for the Gaussian quadrature rule along x direction
num_row = method_data.nquad(2);     % Points for the Gaussian quadrature rule along y direction


% Assemble the rhs which is constant (ext force + Neumann condition)
    b    = op_f_v_tp (sp, msh, f);
    norm_b = sqrt(sum(b.^2));

    % Apply Neumann boundary conditions (ERRORE!!)
    for iside = nmnn_sides
    %Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
      gside = @(varargin) g(varargin{:},iside);
      dofs = sp.boundary(iside).dofs;
     b(dofs) = b(dofs) + op_f_v_tp (sp.boundary(iside), msh.boundary(iside), gside);
    end
    
    %Setto le cond di Dirichlet
    [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp, msh, h, drchlt_sides);
    
    %Initial guess    
    err_d = 1000;
    err_r = 1000;
    num_it=0;
    
    %N-R cycle
    while err_d > eps_d && err_r > eps_r && num_it<num_max_it
             [val, grid] = sp_eval (u, sp, geometry, x, {'value', 'gradient'}); %tilde al posto di grid
              D = val{2};

             f_s  = op_f_d_s_tp(sp, msh, D, num_row, num_col, mat_property);
             rhs  = f_s - b;
             mat    = 0*op_mat_stiff_tp (sp, sp, msh, D, num_row, num_col, mat_property)+op_geo_stiff_tp (sp, sp, msh, D, num_row, num_col, mat_property);
             
             % Apply Dirichlet boundary conditions, for variation they are
             % always 0
             delta_u = zeros (sp.ndof, 1);            
             delta_u(drchlt_dofs) = 0*u_drchlt;
 
             int_dofs = setdiff (1:sp.ndof, drchlt_dofs);
             rhs(int_dofs) = rhs(int_dofs) - mat (int_dofs, drchlt_dofs) * u_drchlt;
 
             % Solve the linearyzed system
             delta_u(int_dofs) = - mat(int_dofs, int_dofs) \ rhs(int_dofs);
            
             norm_delta_u = sqrt(sum(delta_u.^2));
             norm_u = sqrt(sum(u.^2));
             norm_rhs = sqrt(sum(rhs.^2));
             err_d = norm_delta_u/norm_u;
             err_r = norm_rhs/norm_b;
        
             %update the solution
             u = u + delta_u;
             num_it=num_it+1;
             fprintf('It N-R: %d \n', num_it);
             
              vtk_pts = {linspace(0, 1, 20), linspace(0, 1, 20)}; 
              [eu, F] = sp_eval (u, sp, geometry, vtk_pts, {'value'});
              [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

              figure
              quiver (X, Y, squeeze(eu(1,:,:)), squeeze(eu(2,:,:)))
              axis equal tight
              title (['Computed displacement IT: ', num2str(num_it)])
    end
end