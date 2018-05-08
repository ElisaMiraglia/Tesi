% OP_K_MAT: assemble the stiffness matrix K_mat = [K_mat(i,j)].
%
%   mat = op_mat_stiff (spu, spv, msh, D);
%   [rows, cols, values] = op_mat_stiff (spu, spv, msh, D);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   D:     matrix form of the constitutive 4th order tensor (6x6)
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries

% function varargout = op_mat_stiff(spu, spv, msh, D)
% 
%   gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], ...
% 		   msh.nqn, spu.nsh_max, msh.nel);
%   gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], ...
% 		   msh.nqn, spv.nsh_max, msh.nel);
% 
%   ndir = size (gradu, 2);
% 
%   rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
%   cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
%   values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
% 
%   jacdet_weights = msh.jacdet .* msh.quad_weights .* coeff;
%   
%   ncounter = 0;
%   for iel = 1:msh.nel
%     if (all (msh.jacdet(:, iel)))
%       gradu_iel = reshape (gradu(:,:,:,1:spu.nsh(iel),iel), spu.ncomp*ndir, msh.nqn, 1, spu.nsh(iel));
% %       gradu_iel = repmat (gradu_iel, [1,1,spv.nsh(iel),1]);
%       gradv_iel = reshape (gradv(:,:,:,1:spv.nsh(iel),iel), spv.ncomp*ndir, msh.nqn, spv.nsh(iel), 1);
% %       gradv_iel = repmat (gradv_iel, [1,1,1,spu.nsh(iel)]);
% 
%       dim = size(gradu_iel,1);
%       
%       jacdet_iel = reshape (jacdet_weights(:,iel), [1,msh.nqn,1,1]);
%       jacdet_gradu = bsxfun (@times, jacdet_iel, gradu_iel)
%       
%       %First we have to compute B_i and B_j
%       
%       S_jacdet_gradu = sum(bsxfun (@times, jacdet_gradu, S), 1);
%       tmp1 = sum (bsxfun (@times, S_jacdet_gradu, gradv_iel), 1);
%       values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = reshape (sum (tmp1, 2), spv.nsh(iel), spu.nsh(iel));
% 
%       [rows_loc, cols_loc] = ndgrid (spv.connectivity(:,iel), spu.connectivity(:,iel));
%       rows(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = rows_loc;
%       cols(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = cols_loc;
%       ncounter = ncounter + spu.nsh(iel)*spv.nsh(iel);
% 
%     else
%       warning ('geopdes:jacdet_zero_at_quad_node', 'op_K_geo: singular map in element number %d', iel)
%     end
%   end
% 
%   if (nargout == 1 || nargout == 0)
%     varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
%                            values(1:ncounter), spv.ndof, spu.ndof);
%   elseif (nargout == 3)
%     varargout{1} = rows(1:ncounter);
%     varargout{2} = cols(1:ncounter);
%     varargout{3} = values(1:ncounter);
%   else
%     error ('op_K_geo: wrong number of output arguments')
%   end
% 
% end

function mat = op_mat_stiff(spu, spv, msh, F)
  
  F = 2*eye(2)
  mat = spalloc (spv.ndof, spu.ndof, 1);
  
  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradu, 2);

  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
      for idof = 1:spv.nsh(iel)
        ishg = reshape(gradv(:,:,:,idof,iel),spv.ncomp,ndir, []);
        for inode = 1:msh.nqn
           B_i = zeros(3,2,msh.nqn);
           B_i(1:2,1:2,inode)=F.*ishg(:,:,inode);
           B_i(3,:,inode)= [23 23];
        end
        for jdof = 1:spu.nsh(iel)
          jshg = reshape(gradu(:,:,:,jdof,iel),spu.ncomp,ndir, []);
          row_ind = repmat(1:spv.nsh(iel)/2, [1,ndir]);
          col_ind = repmat(1:spu.nsh(iel)/2, [1,ndir]);
          for inode = 1:msh.nqn
             B_j = zeros(3,2,msh.nqn);
             B_j(1:2,1:2,inode)=F.*jshg(:,:,inode);
             B_j(3,:,inode)= [23 23];
             D=eye(3);
             tmp1(:,:,inode) = permute(B_i(:,:,inode), [2 1 3])*D*B_j(:,:,inode);
          end
          
          mat_loc(row_ind(idof), col_ind(jdof)) = mat_loc(row_ind(idof), col_ind(jdof)) + ...
             sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* tmp1(1,1,:));
          
          mat_loc(row_ind(idof), col_ind(jdof)+spu.nsh/2) =  mat_loc(row_ind(idof), col_ind(jdof)+spu.nsh/2) + ...
             sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* tmp1(1,2,:));
          
          mat_loc(row_ind(idof)+spv.nsh/2, col_ind(jdof)) = mat_loc(row_ind(idof)+spv.nsh/2, col_ind(jdof)) + ...
             sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* tmp1(2,1,:));
          
          mat_loc(row_ind(idof)+spv.nsh/2, col_ind(jdof)) = mat_loc(row_ind(idof)+spv.nsh/2, col_ind(jdof)) + ...
             sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* tmp1(2,2,:));
         
        end
      end
      
      mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
        mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) + mat_loc;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradu_gradv: singular map in element number %d', iel)
    end
  end

end



