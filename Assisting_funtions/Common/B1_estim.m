function B1 = B1_estim(radius,ax_dist)

% Constants
mu_0 = 4e-7*pi; % H/m

% Calculations from pag 20 Tobgay 2004
B1 = mu_0*(radius^2)*1/(2*(radius^2+ax_dist^2)^(3/2));

end