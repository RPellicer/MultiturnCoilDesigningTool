function [Bx,By,Bz] = magnetic_field_current_loop_Ruben(x,y,z,x_p,y_p,z_p,a,I0)

%Kilian O'Donoghue
%30th July 2013
% Modified Ruben Pellicer-Guridi CAI, UQ, February 2017
%{
%This function calculates the magnetic field resulting from a single
%circular coil of radius a, carrying a currnet I0. The coil points in the z
%direction. The centre of the coil is located at {x_p,y_p,z_p} and the
%magnetic field is calculated at the point or points stored in the arrays
%x, y and z. This code accepts meshgrid inputs or single values along with
%vectors.

%These equations are based on those given in "General relation for the vector
%magnetic field of a circular current loop: a closer look" by Robert
%Schill, this is available from IEEE.

%% Another way (more accurate) of calculating it with an inhouse equation
%%     rt = sqrt((pointsRef(2,:)+1e-38).^2+pointsRef(3,:).^2); % Centre is along X
%     zt = pointsRef(1,:)+1e-38;
%     at = radi*ones(1,length(zt));
%     %  According to Simson-2001-Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop.pdf
%     alp = sqrt(at.^2+rt.^2+zt.^2-2.*at.*rt);
%     bet = sqrt(at.^2+rt.^2+zt.^2+2.*at.*rt);
%     k3 = sqrt(1-(alp.^2)./(bet.^2)).^2; % IT HAS TO BE SQUARED because "ellipke()" gets k^2 as input intead of k
%     uo = 4*pi*1e-7;          %magnetic permeability (in S.I. units)
%     I = 1;                  %current int the coil (Amp)
%     C = uo*I/pi;
%     [K,E] = ellipke(k3); 
%     Br = (C.*zt./(2.*alp.^2.*bet.*rt)).*((at.^2+rt.^2+zt.^2).*E-(alp.^2).*K); % Eq. 24
%     Bz = (C./(2.*alp.^2.*bet)).*((at.^2-rt.^2-zt.^2).*E+(alp.^2).*K); % Eq. 25
%     [Bx,By] = pol2cart(atan2(pointsRef(2,:),pointsRef(3,:)),Br);

%}
% global u0 %permeability of free space is a global variable
    u0 = 4*pi*1e-7;          %magnetic permeability (in S.I. units) Ruben

    rc=((x-x_p).^2+(y-y_p).^2).^.5; %Radial component is required for cylindrical coordinate system.

    %% Ruben- Equations fail when rc is to small (dibided by a tity number)
    rc_mask = rc>(a/1e6);
    rc=rc.*rc_mask;
    %%
    m=(4.*a.*rc).*(((rc+a).^2)+((z-z_p).^2)).^(-1); %This is a parameter for calculating the Elliptical integrals
    kofkc=(pi/2)+(pi/8).*m+(9*pi/128).*m.^2; %K(k) elliptical function, this is a taylor expansion of the K elliptical integral.
    eofkc=(pi/2)+(-pi/8).*m+(-3*pi/128).*m.^2;%E(k) elliptical function this is a taylor expansion of the E elliptical integral.

    %Note for improved accuracy, Matlab has built in elliptical integral
    %calculation but these expressions here are still very accurate when rc < a

    Brc=(u0.*I0./(2.*pi.*rc)).*(z-z_p).*((((rc+a).^2)+((z-z_p).^2)).^(-.5)).*(-kofkc+eofkc.*((rc.^2+a.^2+(z-z_p).^2)./(((rc-a).^2)+((z-z_p).^2)))); %radial component of B%
    Bz=(u0.*I0./(2.*pi)).*((((rc+a).^2)+((z-z_p).^2)).^(-.5)).*(kofkc-eofkc.*((rc.^2-a.^2+(z-z_p).^2)./(((rc-a).^2)+((z-z_p).^2)))); %axial component of B
    Bx=Brc.*(x-x_p)./rc; %This converts the polar component into cartesian form.
    By=Brc.*(y-y_p)./rc;

    %The following sets any terms that result in Inf to zero, this occurs at
    %the points near the coil itself.
    Bx(isnan(Bx)) = 0 ;
    By(isnan(By)) = 0 ;
    Bz(isnan(Bz)) = 0 ;
end





   


   
