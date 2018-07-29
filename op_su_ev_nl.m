% OP_SU_EV_NL: creates the mat_stiff following linear elasticity path
%
%   mat = op_su_ev_nl (spu, spv, msh, d, num_row, mat_property);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   d:     gradient of the solution one step before
%   mat_prop = vector of 3 elements containing the material properties (Implementation for Mooney Rivlin materials)
%
%
% OUTPUT:
%
%   mat:    assembled matrix

 function mat = op_su_ev_nl (spu, spv, msh, d, num_row, mat_property)
  
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
      
    if (all (msh.jacdet(:,iel)))
      mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
      for idof = 1:spv.nsh(iel)
        ishg  = gradv(:,:,:,idof,iel);
        ishgt = permute (ishg, [2, 1, 3]);
        size(ishgt)
        size(def_grad)
        %ieps  = reshape(ishg + ishgt, spv.ncomp * ndir, [])/2;
        
        ishgtdf = zeros(size(def_grad));
        for inode = 1:size(def_grad, 3)
            ishgtdf(:,:,inode) = ishgt(:,:,inode)*def_grad(:,:,inode);
            ishgdft(:,:,inode) = ishgt(:,:,inode)*def_grad(:,:,inode);
        end
   
        
        for jdof = 1:spu.nsh(iel) 
          jshg  = gradu(:,:,:,jdof,iel);
          jshgt = permute (jshg, [2, 1, 3]);
          jeps  = reshape(jshg + jshgt, spu.ncomp * ndir, [])/2;
          jdiv  = spu.shape_function_divs(:, jdof, iel);
 % The cycle on the quadrature points is vectorized         
          mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
              sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* ...
                   (2 * sum (ieps .* jeps, 1).' .* mu(:,iel)  + ...
                    (idiv .* jdiv) .* lambda(:,iel)));          
        end
      end
      mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
        mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) + mat_loc;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_su_ev: singular map in element number %d', iel)
    end
  end

end
