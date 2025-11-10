function [N_list, Re, tmax, CFL, CK, m, usr] = readInputData()
    % General, numerical and physical
    N_list = [20, 40];
    Re = 40;
    tmax = 1;
    CFL = 0.025;

    % LES
    CK = 0.05;
    m = 2;

    % User Visualization
    usr.displayProgress = true;
    usr.doPlots = false;

end