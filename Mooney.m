%Moonley: assemble the nominal stress S and the matrix form of the constitutive 4th order tensor D of Mooney-Rivlin materials = [K_geo(i,j)].
%
%   [S, D] = Mooney(F, A10, A01, K);
%   
% INPUT:
%
%   F:  Deformation gradient
%   A10, A01: Properties of material
%   K: Bulk modulus (for incompressible material)
%
% OUTPUT:
%
%   S assembled nominal stress matrix
%   D 6x6 matrix that represents the constitutive 4th order tensor C


function varargout = Mooney(F, mat_property)

    A10 = mat_property(1);
    A01 = mat_property(2);
    
    C = F'*F;
    
    %Invariants:
    I1 = sum(diag(C));
    I2 = 1./2 .* ((sum(diag(C))).^2 - sum(diag(C.*C)));
    I3 = det(C);
    
    %Reduced invariants
    %J1= I1 .* I3.^(-1/3);
    %J2= I2.* I3.^(-2/3);
    J3= sqrt(I3);
    
    %Derivatives of I_j invariants of C
    I1E = diag(2.*ones(size(C,1),1));
    I2E = 2.*(diag(I1.*ones(size(C,1),1))-C);
    I3E = 2.* I3 .* C^(-1);
    
    %Derivatives of J_j reduced invariants of C
    J1E = I1E.*(I3)^(-1/3)- 1/3 .* I1.*(I3)^(-4/3).*I3E;
    J2E = I2E.*(I3)^(-2/3)- 2/3 .* I2.*(I3)^(-5/3).*I3E;
    J3E = 1/2 * (I3).^(-1/2).*I3E;
    
    if(length(mat_property)==3)
       K = mat_property(3);
       S = A10.*J1E+A01.*J2E+K.*(J3-1).*J3E;
    else
       S = A10.*J1E+A01.*J2E;
    end 
    
    
  if (nargout == 1)
    varargout{1} = S;
                       
  elseif (nargout == 2)
    varargout{1} = S;
    varargout{2} = D;
  else
    error ('Mooney: wrong number of output arguments')
  end
end