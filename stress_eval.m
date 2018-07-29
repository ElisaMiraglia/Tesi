function sigma_stress = stress_eval (u, sp, geometry, vtk_pts, mat_property)

    [val, ~] = sp_eval (u, sp, geometry, vtk_pts, {'value', 'gradient'});
    u_grad = val{2};
    
    Id = repmat(eye(2), [1,1,size(u_grad,3),size(u_grad,3)]);%Devo creare una matrice identità per numero di nodi volte
    def_grad = bsxfun(@plus, Id, u_grad);
   
    S = zeros(size(u_grad));
    sigma_stress = zeros(size(u_grad));
    for i = 1:size(vtk_pts{1},2)
        for j = 1:size(vtk_pts{2},2) 
           S(:,:,i,j) = Mooney(def_grad(:,:,i,j), mat_property);
           [~,~,sigma_stress(:,:,i,j)] = Mooney(def_grad(:,:,i,j), mat_property);
        end
    end
    
end