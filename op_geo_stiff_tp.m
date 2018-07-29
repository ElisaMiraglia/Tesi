% OP_K_GEO_TP: assemble the stiffness matrix K_geo = [K_geo(i,j)], , exploiting the tensor product structure.
%
%   mat = op_geo_stiff_tp (spu, spv, msh, S);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   d:     gradient of the solution one step before
%   num_col = number of quadrature points in x direction for each element
%   num_row = number of quadrature points in y direction for each element
%   mat_prop = vector of 3 elements containing the material properties (Implementation for Mooney Rivlin materials)
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix

function varargout = op_geo_stiff_tp (space1, space2, msh, d, num_row, num_col, mat_prop, nel_small)

  A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col, 'value', false, 'gradient', true);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, 'gradient', true);
    d_col = d(:,:,(iel-1)*num_col+1:iel*num_col,:);
      
    A = A + op_geo_stiff(sp1_col, sp2_col, msh_col, d_col, num_row, mat_prop, nel_small);
  end

  varargout{1} = A;
  
end
