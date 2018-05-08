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
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% function varargout = op_geo_stiff(spu, spv, msh, S)
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
%       %First we have to compute the mat_vec product 
%       %S calculated outside in order to calculate it only one time.
%       S_d = S(1:dim, 1:dim); %2x2 for bidimensional problem, 3x3 for 3D
%       S_d_rep = repmat(S_d(1,:)', [1, size(jacdet_gradu,2) ]);
%       S__gradu = sum (bsxfun (@times, S_d_rep, jacdet_gradu), 1); 
%       S_d_rep = repmat(S_d(2,:)', [1, size(jacdet_gradu,2) ]);
%       S__gradu(2,:,:,:) = sum (bsxfun (@times, S_d_rep, jacdet_gradu), 1);  
%       
%       if (dim == 3)
%         S_d_rep = repmat(S_d(3,:)', [1, size(jacdet_gradu,2) ]);
%         S__gradu(3,:,:,:) = sum (bsxfun (@times, S_d_rep, jacdet_gradu), 1); 
%       end
%       
%       S_jacdet_gradu = sum(bsxfun (@times, jacdet_gradu, S), 1);
%       tmp1 = sum (bsxfun (@times, S_jacdet_gradu, gradv_iel), 1); %prodotto scalare
%       values(ncounter+(1:spu.nsh(iel)*spv.nsh(iel))) = reshape (sum (tmp1, 2), spv.nsh(iel), spu.nsh(iel)); %Somma nei nodi di quadratura
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

function mat = op_geo_stiff(spu, spv, msh, S)

  mat = spalloc (spv.ndof, spu.ndof, 1);
  
  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], msh.nqn, spu.nsh_max, msh.nel);
  gradv = reshape (spv.shape_function_gradients, spv.ncomp, [], msh.nqn, spv.nsh_max, msh.nel);

  ndir = size (gradu, 2);

  for iel = 1:msh.nel
    if (all (msh.jacdet(:,iel)))
      mat_loc = zeros (spv.nsh(iel), spu.nsh(iel));
      for idof = 1:spv.nsh(iel)
        ishg = reshape(gradv(:,:,:,idof,iel),spv.ncomp, []);
        for jdof = 1:spu.nsh(iel) 
          jshg = reshape(gradu(:,:,:,jdof,iel),spu.ncomp,ndir, []);
          for inode = 1:msh.nqn
             jshg(:,:,inode)=S*jshg(:,:,inode);
           end
             mat_loc(idof, jdof) = mat_loc(idof, jdof) + ...
             sum (msh.jacdet(:,iel) .* msh.quad_weights(:, iel) .* ...
                  sum (ishg .* jshg, 1));  
        end
      end
      mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) = ...
        mat(spv.connectivity(:, iel), spu.connectivity(:, iel)) + mat_loc;
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_gradu_gradv: singular map in element number %d', iel)
    end
  end

end



