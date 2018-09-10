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
%   num_row = number of quadrature points in y direction for each element
%   mat_prop = vector of 3 elements containing the material properties (Implementation for Mooney Rivlin materials)
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix

function mat = op_mat_stiff(spu, spv, msh, d_d, num_row, mat_property, nel_small)
  
  mat = spalloc (spv.ndof, spu.ndof, 1);
  
  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradu, 2);

  for iel = 1:msh.nel
   d_el = d_d(:,:,:,num_row*(iel-1)+1:num_row*iel);
   d_el = reshape(d_el, size(d_el,1), size(d_el,2) ,size(d_el,3)*size(d_el,4));
   
   Id = repmat(eye(2), [1,1,size(d_el,3)]);
  
   def_grad = bsxfun(@plus, Id, d_el);
   
   D = zeros(3,3,size(def_grad,3));
    
   %Mooney Rivlin constitutive law tensor
   for inode = 1:size(def_grad, 3)
%            if (iel<nel_small)
%                 [~,D_node] = Mooney(def_grad(:,:,inode), mat_property);
%            else
%                 [~,D_node] = Mooney(def_grad(:,:,inode), [30*mat_property(1),mat_property(2),20*30*mat_property(1)]);
%            end
           [~,D_node] = Mooney(def_grad(:,:,inode), mat_property);
           D(:,:,inode)=D_node;
   end

%    %%%%LINEAR 
%      E  =  1; nu = .3; 
%     lambda = (nu*E)/((1+nu)*(1-2*nu));
%     mu = E/(2*(1+nu));
%     
%     for inode = 1:size(def_grad, 3)
%                 D_node = [lambda+2*mu lambda 0; lambda lambda+2*mu 0; 0 0 mu];
%                 D(:,:,inode)=D_node;
%     end
%    %%%%%


   row_ind = repmat(1:spv.nsh(iel)/2, [1,ndir]);
   col_ind = repmat(1:spu.nsh(iel)/2, [1,ndir]);
    if (all (msh.jacdet(:,iel)))
      mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
      for idof = 1:spv.nsh(iel)/2
        ishg = reshape(gradv(:,:,:,idof,iel),spv.ncomp,ndir, []);
        B_i = zeros(3,2,msh.nqn);
        for inode = 1:msh.nqn
           ishgsc = sum(ishg(:,:,inode)); %sum over the columns in order to have scalar value functions
           B_i(1:2,1:2,inode)=[def_grad(1,1,inode)*ishgsc(1),def_grad(2,1,inode)*ishgsc(1);def_grad(1,2,inode)*ishgsc(2),def_grad(2,2,inode)*ishgsc(2)];
           B_i(3,:,inode)= [def_grad(1,1,inode)*ishgsc(2)+def_grad(1,2,inode)*ishgsc(1), def_grad(2,1,inode)*ishgsc(2)+def_grad(2,2,inode)*ishgsc(1)];
         end
        for jdof = 1:spu.nsh(iel)/2
          jshg = reshape(gradu(:,:,:,jdof,iel),spu.ncomp,ndir, []);
          B_j = zeros(3,2,msh.nqn);
          tmp1 = zeros(2,2,msh.nqn);
          for inode = 1:msh.nqn
           jshgsc = sum(jshg(:,:,inode)); %sum over the columns in order to have scalar value functions
           B_j(1:2,1:2,inode)=def_grad(:,:,inode).*jshgsc';
           B_j(3,:,inode)= [def_grad(1,1,inode)*jshgsc(2)+def_grad(1,2,inode)*jshgsc(1) def_grad(2,1,inode)*jshgsc(2)+def_grad(2,2,inode)*jshgsc(1)];
           tmp1(:,:,inode) = permute(B_i(:,:,inode), [2 1 3])*D(:,:,inode)*B_j(:,:,inode);
          end
      
          mat_loc(row_ind(idof), col_ind(jdof)) =  mat_loc(row_ind(idof), col_ind(jdof)) + ...
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
    %full(mat)
    %pause
end



