%% this Function helps assign the values to Calculator of Sensitivity (calcSens) 
function s = Final_solver(x,tProbe,visu)
    tProbe.Coil.N_v = x(1);
    tProbe.Coil.N_l = x(2);
    if tProbe.Coil.discreteWire == 1
        tProbe.Coil.do = (x(3) + x(4))*1e-4;    % Ruben 2020 Discretise wire choices
        tProbe.Coil.di = x(4)*1e-4;             % Ruben 2020 Discretise wire choices
    else
        tProbe.Coil.do = x(3) + x(4);    
        tProbe.Coil.di = x(4);
    end
    % Set other needed variables
    tProbe.Coil.Rin = tProbe.Coil.Rout_glob - tProbe.Coil.N_l*tProbe.Coil.do;
    percInsul = 0.99;
    if tProbe.Coil.Rin <= 0
        tProbe.Coil.PenalyzeScaling = tProbe.Coil.PenalyzeScaling*(1+10000*abs(tProbe.Coil.Rin));
        tProbe.Coil.Rin = tProbe.Coil.do;
    end
       
    if (tProbe.Coil.discreteWire == 0) && ((tProbe.Coil.do*percInsul) <= tProbe.Coil.di)
        if tProbe.Coil.do <= tProbe.Coil.di
            tProbe.Coil.PenalyzeScaling = tProbe.Coil.PenalyzeScaling*(1+ 10000*(tProbe.Coil.di - tProbe.Coil.do));
            tProbe.Coil.di = tProbe.Coil.do-0.00001;
        else
            tProbe.Coil.PenalyzeScaling = tProbe.Coil.PenalyzeScaling*(1+ 5*((1/(1-percInsul))*((tProbe.Coil.di-tProbe.Coil.do*percInsul)/tProbe.Coil.di))); % The second part makes the sensitivity 5 times worse when do=di. It can be adjusted
        end
    end
    s = sensitivityI2V(tProbe,visu);
end