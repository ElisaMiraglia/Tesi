%   Moonley: assemble the nominal stress S and the matrix form of the
%   constitutive 4th order tensor D of Mooney-Rivlin materials:
%   Strain energy function expression:
%           W(I1,I2) = A10(I1-3)+A01(I2-3)+K/2*(J-1)^2

%   S = Mooney(F, [A10, A01, K]);
%   [S, D] = Mooney(F, [A10, A01, K]);
%   
%  INPUT:
%
%   F:  Deformation gradient
%   A10, A01: Properties of material
%   K: Bulk modulus (for incompressible material)
%
%  OUTPUT:
%
%   S assembled nominal stress matrix
%   D 6x6 matrix that represents the constitutive 4th order tensor C
%   sigma Cauchy stress tensor


function varargout = Mooney(F, mat_prop)
    
    A10 = mat_prop(1);
    A01 = mat_prop(2);
    K = mat_prop(3);
    
    dim = size(F,1);
    
    if(dim == 2)
        C = eye(3);
        C(1:2, 1:2) = F'*F;
    else
        C = F'*F;
    end
    
    C_i= C^(-1);
    
    C1=C(1,1); C2=C(2,2); C3=C(3,3); C4=C(1,2); C5=C(2,3); C6=C(1,3);
    
    I1 = C1+C2+C3;
    I2 = C1*C2+C2*C3+C1*C3-C4^2-C5^2-C6^2;
    I3 = det(C);
    J3 = sqrt(I3);

    I1E = 2*[1 1 1 0 0 0]';
    I2E = 2*[C2+C3, C3+C1, C1+C2, -C4, -C5, -C6]';
    I3E = 2*[C2*C3-C5^2,  C3*C1-C6^2,  C1*C2-C4^2, ...
         C5*C6-C3*C4, C6*C4-C1*C5, C4*C5-C2*C6]';

   J1E = I3^(-1/3)*I1E - 1/3*I1*I3^(-4/3)*I3E;
   J2E = I3^(-2/3)*I2E - 2/3*I2*I3^(-5/3)*I3E;
   J3E = 1/2*I3^(-1/2)*I3E;
   
   S = A10*J1E + A01*J2E + K*(J3-1)*J3E;

    D=zeros(6);
    I1EE =0;
    I2EE = [0  4  4  0  0  0; 4  0  4  0  0  0; 4  4  0  0  0  0;
         0  0  0 -2  0  0; 0  0  0  0 -2  0; 0  0  0  0  0 -2];
    I3EE = [ 0     4*C3  4*C2  0    -4*C5  0;
          4*C3  0     4*C1  0     0    -4*C6;
          4*C2  4*C1  0    -4*C4  0     0;
          0     0    -4*C4 -2*C3  2*C6  2*C5;
         -4*C5  0     0     2*C6 -2*C1  2*C4;
          0    -4*C6  0     2*C5  2*C4 -2*C2]; 
      
    
    J1EE = -2/3*I3^(-1/2)*(J1E*J3E' + J3E*J1E') + 8/9*I1*I3^(-4/3)*(J3E*J3E') - 1/3*I1*I3^(-4/3)*I3EE;
    J2EE = -4/3*I3^(-1/2)*(J2E*J3E' + J3E*J2E') + 8/9*I2*I3^(-5/3)*(J3E*J3E') + I3^(-2/3)*I2EE - 2/3*I2*I3^(-5/3)*I3EE;
    J3EE = -I3^(-1/2)*(J3E*J3E') + 1/2*I3^(-1/2)*I3EE;


    D = A10*J1EE + A01*J2EE + K*(kron(I3E',I3E)) + K*(I3-1)*J3EE;

     if(dim == 2)
        S_out = [S(1),S(4);S(4),S(2)];
        D_out = zeros(3);
        D_out(1:2,1:3)=D(1:2,[1 2 4]);
        D_out(3,1:3)=D(4,[1 2 4]);
     else
        S_out=S;
        D_out=D;
     end
    
  if (nargout == 1)
    varargout{1} = S_out;
                       
  elseif (nargout == 2)
    varargout{1} = S_out;
    varargout{2} = D_out;
  elseif (nargout ==3)
    varargout{1} = S_out;
    varargout{2} = D_out;
    varargout{3} = 1/det(F)*F*S_out*F';
  else
    error ('Mooney: wrong number of output arguments')
  end
end