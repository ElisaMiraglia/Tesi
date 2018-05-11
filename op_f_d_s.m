% OP_F_D_S: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i).
% 
function rhs = op_f_d_s(spv, msh, d, num_row, mat_property)
 
 rhs   = zeros(spv.ndof, 1);
 gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

 ndir=size(gradv,2);

 for iel = 1:msh.nel
   d_el = d(:,:,:,num_row*(iel-1)+1:num_row*iel);
   d_el = reshape(d_el, size(d_el,1), size(d_el,2) ,size(d_el,3)*size(d_el,4));
   
   Id = repmat(eye(2), [1,1,size(d_el,3)]);%Devo creare una matrice identità per numero di nodi volte
   def_grad = bsxfun(@plus, Id, d_el);
   
   S = zeros(3,size(def_grad,3));
   for inode = 1:size(def_grad, 3)
            S_node = Mooney(def_grad(:,:,inode), mat_property);
            S(:,inode)=[S_node(1,1), S_node(2,2), S_node(1,2)];
   end

   if (all (msh.jacdet(:,iel)))
     rhs_loc = zeros (spv.nsh(iel), 1);
      for idof = 1:spv.nsh(iel)
        row_ind = repmat(1:spv.nsh(iel)/2, [1,ndir]);

        ishg = reshape(gradv(:,:,:,idof,iel),spv.ncomp,ndir, []);
        tmp1=zeros(2,msh.nqn);
        for inode = 1:msh.nqn
           B_i = zeros(3,2,msh.nqn);
           B_i(1:2,1:2,inode)=def_grad(:,:,inode).*ishg(:,:,inode);
           B_i(3,:,inode)= [def_grad(1,1,inode)*sum(ishg(:,2,inode))  def_grad(1,2,inode)*sum(ishg(:,1,inode))];
           tmp1(:,inode) = permute(B_i(:,:,inode), [2 1 3])*S(:,inode);
        end
% The cycle on the quadrature points is vectorized
        %for inode = 1:msh.nqn 
         rhs_loc(row_ind(idof)) = rhs_loc(row_ind(idof)) + ...
           sum(msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .*tmp1(1,:));
         rhs_loc(row_ind(idof)+spv.nsh(iel)/2) = rhs_loc(row_ind(idof)+spv.nsh(iel)/2) + ...
           sum(msh.jacdet(:, iel) .* msh.quad_weights(:, iel) .*tmp1(2,:));
       %end  
     end
     rhs(spv.connectivity(:, iel)) = rhs(spv.connectivity(:, iel)) + rhs_loc; 
   else
     warning ('geopdes:jacdet_zero_at_quad_node', 'op_f_d_s: singular map in element number %d', iel)
   end
 end
 
end
