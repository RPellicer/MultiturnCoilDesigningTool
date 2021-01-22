function B1 = B1_estimMultiTurn(pos,ax_dist)
% "pos" should contain the position from the centre of the turn to the axis and the elevation of each loop. 
% "ax_dist" is the distance from the surface of the coil to the ROI

% Constants
mu_0 = 4e-7*pi; % H/m

radius = squeeze(pos(:,1));
ax_dist = ax_dist + squeeze(pos(:,2));

% Calculations from pag 20 Tobgay 2004
B1M = mu_0*(radius.^2)*1./(2*(radius.^2+ax_dist.^2).^(3/2));
B1 = sum(sum(B1M));
end