% OP_F_D_S_TP: f_d(i) = int(B_i^T, S(d)), exploiting the tensor product structure.
%
%   S(d) = (S11, S22, S12)
%   f = op_f_d_s_tp (spv, msh, d);
%
% INPUT:
%     
%   spv:   object representing the function space (see sp_vector)
%   msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%   d:     gradient of the solution one step before
%   num_col = number of quadrature points in x direction for each element
%   num_row = number of quadrature points in y direction for each element
%   mat_prop = vector of 3 elements containing the material properties (Implementation for Mooney Rivlin materials)

%
% OUTPUT:
%
%   f: assembled part of right-hand side of linearized problem

function rhs = op_f_d_s_tp (space, msh, d, num_col, num_row, mat_prop,nel_small)
  rhs = zeros (space.ndof, 1);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col = sp_evaluate_col (space, msh_col, 'value', false, 'gradient', true);
    d_col = d(:,:,(iel-1)*num_col+1:iel*num_col,:);

    rhs = rhs + op_f_d_s(sp_col, msh_col, d_col, num_row, mat_prop, nel_small);
  end

end
