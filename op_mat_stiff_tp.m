% OP_K_MAT_TP: assemble the stiffness matrix K_mat = [K_mat(i,j)], , exploiting the tensor product structure.
%
%   mat = op_geo_mat_tp (spu, spv, msh, D);
%   [rows, cols, values] = op_geo_stiff_tp (spu, spv, msh, S); (DA IMPLEMENTARE)
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   D:     matrix form of the constitutive 4th order tensor (6x6 for 3D problems, 3x3 for 2D problems)
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries (da implementare)
%   cols:   column indices of the nonzero entries (da implementare)
%   values: values of the nonzero entries (da implementare)

function varargout = op_mat_stiff_tp (space1, space2, msh, D)

  A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp1_col = sp_evaluate_col (space1, msh_col, 'value', false, 'gradient', true);
    sp2_col = sp_evaluate_col (space2, msh_col, 'value', false, 'gradient', true);

    for idim = 1:msh.rdim
         x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end
    coeffs = D(x{:});

  A = A + op_geo_stiff(sp1_col, sp2_col, msh_col, coeffs);
  end

  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end
end