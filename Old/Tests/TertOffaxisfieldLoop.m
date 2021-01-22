clear all,  clc
axiallimit = 2.0; % meters from center
radiallimit = 0.5; % maximum radius to investigate
curveqty = 5;
z = linspace(0,axiallimit);
rvec = linspace(0, radiallimit, curveqty);

a =1;   % Coil radius
% z =3;   % Distance from the coil along the axis 
for i = 1 : length(rvec)
    r =rvec(i);   % Off axis distance

    I =1;   % Current
    uo = 4*pi*1e-7;
%%############################################ Encoded
    at =a;   % Coil radius
    zt =z;   % Distance from the coil along the axis 
    rt =r;   % Off axis distance

    k2 = abs(4*rt.*at ./(zt.^2 +(at+rt).^2));        %input for elliptic integrals
    [K,E] = ellipke(k2);           %elliptic integrals of first and second kind 
    Bz2 = (((uo*I)/(2*pi))*(1./sqrt(zt.^2 +(at+rt).^2)).*((((at.^2) - (zt.^2) - (rt.^2))./(zt.^2 + (rt-at).^2)).*E + K));
%%############################################ http://nbviewer.jupyter.org/github/tiggerntatie/emagnet-py/blob/master/offaxis/off_axis_loop.ipynb
    x = z;
    % Loop magnetic field calculus by using Biot-Savart law with respect elliptic integrals.pdf
    B0 = uo*I/(2*a);
    alpha = r./a;
    betha = x./a;
    gamma = x./r;
    Q = (1+alpha).^2+betha.^2; 
    k = sqrt(4.*alpha./Q);
    [K,E] = ellipke(k); 
    Bz1 = B0./(pi.*sqrt(Q)).*(E.*(1-alpha.^2-betha.^2./(Q-4.*alpha))+K);
%%############################################ "Simple Analytic Expressions
%%for the Magnetic Field of a Circular Current Loop.pdf"
alp = sqrt(a.^2+r.^2+z.^2-2*a*r);
bet = sqrt(a.^2+r.^2+z.^2+2*a*r);
k3 = sqrt(1-(alp.^2)./(bet.^2)).^2;
C = uo*I/pi;
    [K,E] = ellipke(k3); 
Bz3 = C./((2.*alp.^2).*bet).*((a.^2-r.^2-z.^2).*E+(alp.^2).*K);
%%############################################ 
    figure(2)
    plot(z,Bz2)
    hold on
    %     plot(z,Bz1,'r')
    %     hold off
    grid on
end
% legend('NewWeb','AlreadyCoded')