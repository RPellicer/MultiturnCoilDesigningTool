function [R_ac_prox,R_ac]=CoilResistACLitz_Commercial(R_dc,rc,zc,di_strand,s,freq,NS,do)
    %% Calculation of the AC resistance
    % From "Savukov -Detection of NMR signals with a radio-frequency atomic magnetometer"
    % Also see "https://www.newenglandwire.com/products/litz-wire-and-formed-cables/theory"
    %% Definition of material properties
    % Copper
    mu_c = 4*pi*1e-7;% absolute magnetic permeability of the conductor [H/m] of an air core inductor
    ro_c = 1.72e-8; % resistivity of the wire

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
    R_ac_prox = R_dc*(1+F+G_2);

    % Comparing it with https://www.newenglandwire.com/traditional-litz-wire-theory/

    % X = [0, 0.5,    0.6,    0.7,    0.8,    0.9,    1     ];
    % H = [1, 1.0003, 1.0007, 1.0012, 1.0021, 1.0034, 1.005];
    % X = 0.2178678*di_strand*2,54e-5*sqrt(freqs) %  CAUTION DI in the equations
    % are in 'mils' (1/1000th of an inch) so they needed to be scaled to meters by 1/2,54e-5

    Htmp = 1; % this is a simplification as the range of H varies only by a 0.5%
    if NS<=3
        Ktmp=1.55;
    elseif NS<=9
        Ktmp=1.84;
    elseif NS<=27
        Ktmp =1.92;
    else
        Ktmp=2;    
    end
    Gtmp = (di_strand*1e3/25.4*sqrt(freq)./10.44).^4;
    RAC2RDC_ratio = Htmp + Ktmp*(NS*di_strand/do)^2.*Gtmp;
    R_ac = R_dc*RAC2RDC_ratio;
end