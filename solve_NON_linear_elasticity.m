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

% Assemble the rhs which is constant (ext force + Neumann condition)
    b    = op_f_v_tp (sp, msh, f);

    % Apply Neumann boundary conditions
    for iside = nmnn_sides
    % Restrict the function handle to the specified side, in any dimension, gside = @(x,y) g(x,y,iside)
      gside = @(varargin) g(varargin{:},iside);
      dofs = sp.boundary(iside).dofs;
      b(dofs) = b(dofs) + op_f_v_tp (sp.boundary(iside), msh.boundary(iside), gside);
    end
    
    symm_dofs = [];
    
    %Initial guess
    u = zeros (sp.ndof, 1);
   
    %N-R cycle
    while(err_d > eps_d && err_r > eps_r)
        
     [grad_u, F] = sp_eval (u, space, geometry, vtk_pts, {'gradient'});
      
     %grad_u è valuatato sulla griglia
        def_grad = eye(2)+grad_u{1}(1,1,:,:);
        [S,D] = Mooney(def_gradient, mat_prop); %MODIFICARE MOONEY PER INODE
        f_s  = op_f_d_s_tp(sp, sp, msh, S); %SCRIVERE LA FUNZIONE
        rhs  = f_s - b;
        mat    = op_geo_stiff_tp (sp, sp, msh, S) + op_mat_stiff_tp (sp, sp, msh, D); 
        
        % Apply Dirichlet boundary conditions
        delta_u = zeros (sp.ndof, 1);  
        [u_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (sp, msh, h, drchlt_sides);
        delta_u(drchlt_dofs) = u_drchlt;

        int_dofs = setdiff (1:sp.ndof, [drchlt_dofs, symm_dofs]);
        rhs(int_dofs) = rhs(int_dofs) - mat (int_dofs, drchlt_dofs) * u_drchlt;

        % Solve the linear system
        delta_u(int_dofs) = mat(int_dofs, int_dofs) \ rhs(int_dofs);
        
        u_new = u + delta_u;
        u = u_new;
    end