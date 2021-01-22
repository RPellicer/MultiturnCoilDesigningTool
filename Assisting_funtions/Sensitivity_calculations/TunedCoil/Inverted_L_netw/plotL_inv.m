%% This function plots the bode of a 'sym' function along the 's' variable
function [h] = plotL_inv(f,vars,vals)
    f_sol = subs(f,vars,vals); 
    [nsym,dsym]=numden(f_sol);
    % convert them to MATLAB polynomials
    mySys=tf(sym2poly(nsym),sym2poly(dsym));
    % Bode Plot
    opts = bodeoptions;
    opts.FreqUnits = 'Hz';
    opts.MagUnits = 'abs';
    bodefreq = logspace(log10(2*pi*9.7e3),log10(2*pi*1.03e4),10000);
    h = bodeplot(mySys,bodefreq,opts);
    grid on
end