fprintf('Adding to the path GeoPDES folders\n')
addpath(genpath('C:/Program Files/MATLAB/R2017b/geopdes'))
addpath(genpath('C:/Program Files/MATLAB/R2017b/nurbs'))
addpath(genpath('C:/Users/utente/desktop/TESI/NonLinear Elasticity'))
fprintf('Done \n')

choice = input('Do you want to go to examples directory? Y/N ', 's');
if(choice=='Y' || choice == 'y')
    cd('C:/Program Files/MATLAB/R2017b/geopdes/inst/examples')
end