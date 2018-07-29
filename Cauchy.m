%   Moonley: assemble the nominal stress S and the matrix form of the
%   constitutive 4th order tensor D of Mooney-Rivlin materials:
%   Strain energy function expression:
%           W(I1,I2) = A10(I1-3)+A01(I2-3)

%   sigma = Mooney(F, [A10, A01, K])
%
%   
%  INPUT:
%
%   F:  Deformation gradient
%   A10, A01: Properties of material
%   K: Bulk modulus (for incompressible material)
%
%  OUTPUT:
%
%   sigma: Cauchy stress tensor

function Cauchy(F, mat_prop)
    
    A10 = mat_prop(1);
    A01 = mat_prop(2);
    K = mat_prop(3);
    
    J =det(F);
    B = F*F';
    I1 = trace(B);
    I2 = 1/2*((trace(B))^2-trace(B^2));
    %p_star = 2/3*(A10*I1-A01*I2);
    %sigma = -p_star.*eye(3)+2*A10.*B-2*A01.*B^(-1);

    I1bar = J^(-2/3)*I1;
    I2bar=J^(-4/3)*I2;
    Bbar=det(B)^(-1/3)*B;
    p=-K*(J-1);
    sigma = 1/J*(-p*eye(3)+2*(A10+I1bar*A01)*Bbar-2*A01*Bbar*Bbar-2/3*(A10*I1bar+2*A01*I2bar)*eye(3));
    
    
    S= Mooney(F,mat_prop);
    
    cauchy_mio = 1/J*F*S*F';
    
    sigma-cauchy_mio
    S
    S_wiki=F^(-1)*J*sigma*F'^(-1)
    
  end
