%% emf for the near field ROI
function [emf_near] = emfNear(freqs,rCoil,hCoil,r_ROI,hBase_ROI,hfinal_ROI,Bp)
    w = 2*pi*freqs;
    VolumePhant = abs(hfinal_ROI-hBase_ROI)*pi*r_ROI^2;
    sens_vol = Bdet_ROI_averagV2(hCoil,rCoil,r_ROI,hBase_ROI,hfinal_ROI);
%%     emf = d/dt int_vol(Bdet(r)*M(r,t)*dvol), from Lang, 'Principles of MRI' pag 95
% or Hoult 1976 'SIGNAL-TO-NOISE OF NUCLEAR RESONANCE EXPERIMENT'  pag 3 (73 of the .pdf)
    emf_near = 1j*w*sens_vol*netMagnetization(Bp)*VolumePhant;
%% ALTERNATIVE Derived from 'Calculated signal-to-noise ratio of MRI detected with SQUIDs and Faraday detectors in fields from 10 lT to 1.5 T'
%     emf_near = 1j*w*sens_vol*netMagnetization(Bp)*VolumePhant;

end

% Test code
% emfNear(freqs,rCoil,hCoil,r_ROI,hBase_ROI,hfinal_ROI,Bp)