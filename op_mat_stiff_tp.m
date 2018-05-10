% OP_K_MAT_TP: assemble the stiffness matrix K_mat = [K_mat(i,j)], , exploiting the tensor product structure.
%
%   mat = op_mat_stiff_tp(space1, space2, msh, d, num_row, num_col, mat_prop)
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   
% OUTPUT:
%
%   mat:    assembled stiffness matrix

function A = op_mat_stiff_tp (space1, space2, msh, d, num_row, num_col, mat_prop)

  A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);
  
  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col, 'value', false, 'gradient', true);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, 'gradient', true);
    d_col = d(:,:,(iel-1)*num_col+1:iel*num_col,:);
    
    A = A + op_mat_stiff(sp1_col, sp2_col, msh_col, d_col, num_row, mat_prop);
  end
end