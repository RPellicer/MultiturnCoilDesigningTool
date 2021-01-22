%% The fuction for average magnetic field (sensitivity) over Region Of Interest (ROI)
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

% This is the flux accros the coil (remaining to multiply by the magnetization) when estimating the EMF!!
function [sens] = Bdet_ROI_averag(hCoil,rCoil,r_ROI,hbase_ROI,hfinal_ROI)
    if length(hCoil)~=length(rCoil)
        sens = 'error: Ambiguity in number of loops';
        return    
    end
    uo = 4*pi*1e-7;          %magnetic permeability (in S.I. units)
    I = 1;                  %current int the coil (Amp)
    r = linspace(0,r_ROI,10);        %span of radius RING STEP
    if hbase_ROI~=hfinal_ROI
        heigh_ROI = linspace(hbase_ROI,hfinal_ROI,10); % 100 HEIGHT STEP
    else
        heigh_ROI=hfinal_ROI;
    end
    B_ROI_disc = zeros(1,length(heigh_ROI));

    for idx=1:length(heigh_ROI)       % Splitting the total height into 100 disks

        z = heigh_ROI(idx) + hCoil;      %normalize the height for the disk considered
        %  According to Simson-2001-Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop.pdf
        rt = repmat(r,length(rCoil),1);
        at = repmat(rCoil',1,length(r));
        zt = repmat(z',1,length(r));
        alp = sqrt(at.^2+rt.^2+zt.^2-2.*at.*rt);
        bet = sqrt(at.^2+rt.^2+zt.^2+2.*at.*rt);
        k3 = sqrt(1-(alp.^2)./(bet.^2)).^2; % IT HAS TO BE SQUARED because "ellipke()" gets k^2 as input intead of k
        C = uo*I/pi;
% To be optimised        
%         m=(4.*a.*rc).*(((rc+a).^2)+((z-z_p).^2)).^(-1); %This is a parameter for calculating the Elliptical integrals
%         kofkc=(pi/2)+(pi/8).*m+(9*pi/128).*m.^2; %K(k) elliptical function, this is a taylor expansion of the K elliptical integral.
%         eofkc=(pi/2)+(-pi/8).*m+(-3*pi/128).*m.^2;%E(k) elliptical function this is a taylor expansion of the E elliptical integral.     
        
        [K,E] = ellipke(k3); 
        Bz = C./((2.*alp.^2).*bet).*((at.^2-rt.^2-zt.^2).*E+(alp.^2).*K);
        Bz_buf(idx,:,:) = Bz;
%         figure(2)
        Bz_circ_tot = sum(Bz,1);     % axial magnetization is summed for all the coils
        % calculation of average magnetic field from the disc
        B_ROI_disc(idx) = 2*pi*trapz(r,(Bz_circ_tot.*r)')/(pi*r_ROI^2); % Average_Bz = doubleCylindricalIntegral(Bz)/Disc_Area ; 
                                                        % doubleCylindricalIntegral(Bz)= 2*pi*int(Bz*radious)
    end
    sens = mean(B_ROI_disc);
%     h = figure;
%     plot(heigh_ROI,BzDist,'-*');
%     title('Axial component of unit coil (1m radius, 1A current) B field for various measurement radii');
%     grid on
end
% To test the correctness of this code uncomment the plot and run the next
% command:
% >> Bdet_ROI_averag(0,1,0.5,0,2) 
% Test function call should look like the plot in this page:
% http://nbviewer.jupyter.org/github/tiggerntatie/emagnet-py/blob/master/offaxis/off_axis_loop.ipynb
