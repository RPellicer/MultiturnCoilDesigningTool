function [R_ac,F,G_2]=CoilResistACLitz(R_dc,rc,zc,di_strand,s,freq,NS)
%% Calculation of the AC resistance
% From "Savukov -Detection of NMR signals with a radio-frequency atomic magnetometer"
% Also see "https://www.newenglandwire.com/products/litz-wire-and-formed-cables/theory"
%% Definition of material properties
% Copper
mu_c = 4*pi*1e-7;% absolute magnetic permeability of the conductor [H/m] of an air core inductor
ro_c = 1.7e-8; % resistivity of the wire

omega = 2*pi*freq;
% s = 2*a; % s is the spacing between the centers of the wires
rc_s = rc/s; % because u is calculated in units of s.
zc_s = zc/s; % because u is calculated in units of s.
rc_s = rc_s-min(rc_s)+1;
zc_s = zc_s-min(zc_s)+1;
rc_s = round(rc_s);
zc_s = round(zc_s);
N_t = length(rc);
% skin depth
skinDepth = sqrt(2*ro_c./(omega.*mu_c));
z = di_strand./(skinDepth*sqrt(2));

% Skin depth EFFECT
F = -(z.^2./8).*imag(besselj(3,z*sqrt(-1i))./besselj(1,z*sqrt(-1i)));

% Proximity EFFECT
G = -(z.^2./8).*imag(besselj(2,z*sqrt(-1i))./besselj(0,z*sqrt(-1i)));
temp_ul = ones(N_t,1)*rc_s;
temp_uv = ones(N_t,1)*zc_s;
mask = ~diag(ones(N_t,1));
Au = (temp_ul - temp_ul')./((temp_ul - temp_ul').^2+(temp_uv - temp_uv').^2);
Bu = (temp_uv - temp_uv')./((temp_ul - temp_ul').^2+(temp_uv - temp_uv').^2);
u = (1/N_t)*(nansum(nansum(mask.*Au).^2)+...
    nansum(nansum(mask.*Bu).^2));

G_2 = (u+2)*(NS*di_strand/s).^2*G;

% Total AC resistance
R_ac = R_dc*(1+F+G_2);
end