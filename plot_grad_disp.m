function [] = plot_grad_disp(eu,F, sigma_stress, S)

[X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));

figure
subplot(1,2,1)
quiver (X, Y, squeeze(eu{1}(1,:,:)), squeeze(eu{1}(2,:,:)))
title ('Displacement')

subplot(1,2,2)
surf(X,Y,sqrt((squeeze(eu{1}(1,:,:))).^2+squeeze(eu{1}(2,:,:).^2)))
view(2); shading interp; colorbar;
title ('Norm')


figure
subplot(2,2,1)
surf (X, Y, squeeze(eu{2}(1,1,:,:)))
view(2); colorbar;
title ('Gradient component \nabla_{xx}')

subplot(2,2,2)
surf (X, Y, squeeze(eu{2}(1,2,:,:)))
view(2); shading interp; colorbar;
title ('Gradient component \nabla_{xy}')

subplot(2,2,3)
surf (X, Y, squeeze(eu{2}(2,1,:,:)))
view(2); shading interp; colorbar;
title ('Gradient component \nabla_{yx}')

subplot(2,2,4)
surf (X, Y, squeeze(eu{2}(2,2,:,:)))
view(2); shading interp; colorbar;
title ('Gradient component \nabla_{yy}')

figure
subplot(2,2,1)
surf (X, Y, squeeze(S(1,1,:,:)))
view(2); shading interp; colorbar;
title ('2-nd Piola Kirchhoff stress component S_{xx}')

subplot(2,2,2)
surf (X, Y, squeeze(S(1,2,:,:)))
view(2); shading interp; colorbar;
title ('2-nd Piola Kirchhoff stress component S_{xy}')


subplot(2,2,3)
surf (X, Y, squeeze(S(2,1,:,:)))
view(2); shading interp; colorbar;
title ('2-nd Piola Kirchhoff stress component S_{yx}')


subplot(2,2,4)
surf (X, Y, squeeze(S(2,2,:,:)))
view(2); shading interp; colorbar;
title ('2-nd Piola Kirchhoff stress component S_{yy}')

end