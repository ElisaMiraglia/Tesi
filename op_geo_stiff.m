% OP_K_GEO: assemble the stiffness matrix K_geo = [K_geo(i,j)].
%
%   mat = op_geo_stiff (spu, spv, msh, S);
%   [rows, cols, values] = op_geo_stiff (spu, spv, msh, S);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_scalar/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_scalar/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col)
%   S:     nominal stress matrix (2x2 for 2D problems, 3x3 for 3D)
%
% OUTPUT:
%
%   mat:    assembled stiffness matrix

function mat = op_geo_stiff(spu, spv, msh, d, num_row, mat_property)

  mat = spalloc (spv.ndof, spu.ndof, 1);
  
  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradu, 2);
  
  for iel = 1:msh.nel
    d_el = d(:,:,:,num_row*(iel-1)+1:num_row*iel);
    d_el = reshape(d_el, size(d_el,1), size(d_el,2) ,size(d_el,3)*size(d_el,4));
   
    Id = repmat(eye(2), [1,1,size(d_el,3)]);
    def_grad = bsxfun(@plus, Id, d_el);
   
    S = zeros(2,2,size(def_grad,3));
        for inode = 1:size(def_grad, 3)
            S_node = Mooney(def_grad(:,:,inode), mat_property);
            S(:,:,inode)=S_node;
        end
   
    if (all (msh.jacdet(:,iel)))
      mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
      for idof = 1:spv.nsh(iel)
        ishg = reshape(gradv(:,:,:,idof,iel),spv.ncomp,ndir, []);
        for jdof = 1:spu.nsh(iel) 
        jshg = reshape(gradu(:,:,:,jdof,iel),spu.ncomp,ndir, []);
        tmp1 = zeros(msh.nqn,1);
          for inode = 1:msh.nqn
             jshg(:,:,inode)=S(:,:,inode)*jshg(:,:,inode);
             tmp1(inode)=sum(sum(ishg(:,:,inode).*jshg(:,:,inode),1));
          end
        mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
             sum (msh.jacdet(:,iel) .* ...
                  tmp1 .* msh.quad_weights(:, iel));
        end
      end
      mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
        mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) + mat_loc;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradu_gradv: singular map in element number %d', iel)
    end
  end

