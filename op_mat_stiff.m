% OP_K_MAT: assemble the stiffness matrix K_mat = [K_mat(i,j)].
%
%   mat = op_mat_stiff(spu, spv, msh, d, num_row, mat_property)
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

function mat = op_mat_stiff(spu, spv, msh, d, num_row, mat_property)
  
  mat = spalloc (spv.ndof, spu.ndof, 1);
  
  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradu, 2);

  for iel = 1:msh.nel
   d_el = d(:,:,:,num_row*(iel-1)+1:num_row*iel);
   d_el = reshape(d_el, size(d_el,1), size(d_el,2) ,size(d_el,3)*size(d_el,4));
   
   Id = repmat(eye(2), [1,1,size(d_el,3)]);%Devo creare una matrice identità per numero di nodi volte
   def_grad = bsxfun(@plus, Id, d_el);
   
   D = zeros(3,3,size(def_grad,3));
   for inode = 1:size(def_grad, 3)
            [S_node,D_node] = Mooney(def_grad(:,:,inode), mat_property);
            D(:,:,inode)=D_node;
   end
   
   row_ind = repmat(1:spv.nsh(iel)/2, [1,ndir]);
   col_ind = repmat(1:spu.nsh(iel)/2, [1,ndir]);
    if (all (msh.jacdet(:,iel)))
      mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
      for idof = 1:spv.nsh(iel)
        ishg = reshape(gradv(:,:,:,idof,iel),spv.ncomp,ndir, []);
        for inode = 1:msh.nqn
           B_i = zeros(3,2,msh.nqn);
           B_i(1:2,1:2,inode)=def_grad(:,:,inode).*ishg(:,:,inode);
           B_i(3,:,inode)= [23 23];
        end
        for jdof = 1:spu.nsh(iel)
          B_j = zeros(3,2,msh.nqn);
          jshg = reshape(gradu(:,:,:,jdof,iel),spu.ncomp,ndir, []);
          for inode = 1:msh.nqn
             B_j(1:2,1:2,inode)=def_grad(:,:,inode).*jshg(:,:,inode);
             B_j(3,:,inode)= [23 23];
             tmp1(:,:,inode) = permute(B_i(:,:,inode), [2 1 3])*D(:,:,inode)*B_j(:,:,inode);
          end
          
          mat_loc(row_ind(idof), col_ind(jdof)) = mat_loc(row_ind(idof), col_ind(jdof)) + ...
             sum(msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* reshape(tmp1(1,1,:), msh.nqn, 1));
         
          mat_loc(row_ind(idof), col_ind(jdof)+spu.nsh/2) =  mat_loc(row_ind(idof), col_ind(jdof)+spu.nsh/2) + ...
             sum(msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* reshape(tmp1(1,2,:), msh.nqn, 1));
         
          mat_loc(row_ind(idof)+spv.nsh/2, col_ind(jdof)) = mat_loc(row_ind(idof)+spv.nsh/2, col_ind(jdof)) + ...
             sum(msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* reshape(tmp1(2,1,:), msh.nqn, 1));
          
          mat_loc(row_ind(idof)+spv.nsh/2, col_ind(jdof)+spu.nsh/2) = mat_loc(row_ind(idof)+spv.nsh/2, col_ind(jdof)+spu.nsh/2) + ...
             sum(msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* reshape(tmp1(2,2,:), msh.nqn, 1));
         
        end
      end
      
      mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
        mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) + mat_loc;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_mat_stiff: singular map in element number %d', iel)
    end
  end

end



