clear all
n = 1000;
err_diag = [];
err_xz = [];
err_yz = [];
err_diag_mio=[];
err_yz_mio = [];
err_xz_mio = [];

for j = 1:n
    F = eye(3);
    F(1:2, 1:2) = rand(2);
    [S2, D2] = Mooney_Online(F, 20, 80, 100);
    [S1, D1] = Mooney(F, [20, 80, 100]);
    
    err_diag = [err_diag, S2(3)];
    err_diag_mio = [err_diag, S1(3,3)];
    
    err_yz = [err_yz, S2(5)];
    err_yz_mio = [err_yz_mio, S1(2,3)];
    
    err_xz = [err_xz, S2(6)];
    err_xz_mio = [err_xz_mio, S1(1,3)];
    
    j = j+1;
end

figure(1)
subplot(1,3,1)
boxplot(err_diag)
title('Componente S33');
subplot(1,3,2)
boxplot(err_yz)
title('Componente S23');
subplot(1,3,3)
boxplot(err_xz)
title('Componente S13');

figure(2)
subplot(1,3,1)
boxplot(err_diag_mio)
title('Errore S33 mio');
subplot(1,3,2)
boxplot(err_yz_mio)
title('Errore S23 mio');
subplot(1,3,3)
boxplot(err_xz_mio)
title('Errore S13 mio');
