%% The fuction for Total magnetization (sensitivity) over Region Of Interest (ROI)
% 18th February,2016, Centre for Advanced Imaging, UQ
% %% INPUT &&------------------------------------------------
% hL        the vector of the position of loops
% A         the vector of the radius wrt to hL
% R         the radius of the disk (extreme point) considered
% hbase     the height of the base disk
% hf         the total height of the disk
% %% OUTPUT %%-----------------------------------------------
% vol       the net volume over the height of the considered space
%           which is net avg magnetization
%--------------------------------------------------------------------------
function [sens] = ROI_sens_total(hL,A,R,hbase,hf)
if length(hL)~=length(A)
    sens = 'error: Ambiguity in number of loops';
    return    
end
uo = 4*pi*1e-7;          %magnetic permeability (in S.I. units)
I = 1;                  %current int the coil (Amp)
r = 0:(R/100):R;        %span of radius

if hbase~=hf
    step = (hf - hbase)/100;
    hv=hbase:step:hf;
else
    hv=hf;
end
V = zeros(1,length(hv));

for j=1:length(hv)       % Splitting the total height into 100 disks
%     Brp = zeros(1,101);
%     for k=1:length(A)      % Loop for each value of Z wrt to a 
%         z = hv(j) + hL(k);      %normalize the height for the disk considered
%         a = A(k);
%     for k=1:length(A)      % Loop for each value of Z wrt to a 
        z = hv(j) + hL;      %normalize the height for the disk considered
        a = A;
%         Bz_circ_tot = zeros(length(a),length(r));
        % Magnetization
%         for i =1:length(r)       %splitting the R in 100 sub-parts
%             k2 = abs(4*r(i)*a /(z^2 +(a+r(i))^2));         %input for elliptic integrals
%             [K,E] = ellipke(k2);           %elliptic integrals of first and second kind 
%             Bz = ((u*I)/(2*pi))*(1/sqrt(z^2 +(a+r(i))^2))*((((a^2) - (z^2) - (r(i)^2))/(z^2 + (r(i)-a)^2))*E + K);
%         %    Br = ((u*I*z)/(2*pi*r))*(1/(sqrt(z^2 +(a+r)^2)))*(((a^2 + z^2 + r^2)/(z^2 + (r-a)^2))*E - K);
%             Bz_circ = Bz*r(i);
%             Bz_circ_tot(1,i) = Bz_circ_tot(1,i)+ Bz_circ;           %axial magnetization
%         %    Brp(1,c) = Brp(1,c)+ Br;         %radial magnetization (use if required)
%          end
            rt = repmat(r,length(a),1);
            at = repmat(a',1,length(r));
            zt = repmat(z',1,length(r));
            k2 = abs(4*rt.*at ./(zt.^2 +(at+rt).^2));         %input for elliptic integrals
            [K,E] = ellipke(k2);           %elliptic integrals of first and second kind 
            Bz = ((uo*I)/(2*pi))*(1./sqrt(zt.^2 +(at+rt).^2)).*((((at.^2) - (zt.^2) - (rt.^2))./(zt.^2 + (rt-at).^2)).*E + K);
        %    Br = ((u*I*z)/(2*pi*r))*(1/(sqrt(z^2 +(a+r)^2)))*(((a^2 + z^2 + r^2)/(z^2 + (r-a)^2))*E - K);
            Bz_circ_tot = sum(Bz,1);     % axial magnetization is added per coil
        %    Brp(1,c) = Brp(1,c)+ Br;         %radial magnetization (use if required)
%     end
    % calculation of average magnetic field from the disc 
    V(j) = 2*pi*trapz(r,(Bz_circ_tot.*rt)')/(pi*R^2); % Average_Bz = doubleCylindricalIntegral(Bz)/Disc_Area ; 
                                                    % doubleCylindricalIntegral(Bz)= 2*pi*int(Bz*radious)
%         % calculation of TOTAL magnetic field from the disc 
%     V(j) = 2*pi*trapz(r,(Bz_circ_tot.*rt)');

end
% p = hbase:((h-hbase)/100):h;
% vol = trapz(p,V);
sens = mean(V);
% sens = sum(V); % Sum from all the discs
end
% you can make few changes concerning the arrays of the height and radius
% of the current carrying loop
% The current in the loop is considered 1 amps for the simple calculations
% of sensitivity