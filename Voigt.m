function [S_out, D_out] = Voigt(S, D, dim)
    S_out = zeros(dim);   
    if dim==3
        S_out = S_out + diag([S(1), S(2), S(3)]);
        S_out(2,3) = S(5);
        S_out(1,3) = S(6);
        S_out(1,2) = S(4);
        S_out(3,1) = S(6);
        S_out(2,1) = S(4);
        S_out(3,2) = S(5);
         
        D_out = D;
    else
        S_out = S_out + diag([S(1), S(2)]);
        S_out(1,2) = S(4);
        S_out(2,1) = S(4);
        
        D_out = zeros(3);
        D_out(1,1) = D(1,1);
        D_out(2,2) = D(2,2);
        D_out(3,3) = D(6,6);
        D_out(1,2) = D(1,2);
        D_out(2,1) = D(1,2);
        D_out(1,3) = D(1,6);
        D_out(3,1) = D(1,6);
        D_out(2,3) = D(2,6);
        D_out(3,2) = D(2,6);
    end
end