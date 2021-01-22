% compare coils
round1 = load('optimResults-RmatchDiscreteFix_Miguel-Janesek8x-50X50mmCylinder-75k_5kBw_v1.mat');
% round2 = load('optimResults-Miguel-Janesek8x-50X50mmCylinder-75k_v3.mat');

tmp = round1.ProbeOpt(5).Sens.SensMean
tmp = round1.ProbeOpt(5).Sens.Sens_peak
tmp = round1.ProbeOpt(5).Sens.SensScore
tmp = round1.ProbeOpt(5).Sens.PenalyzeScaling
tmp = round1.ProbeOpt(5).Coil.SelfRes

[s_test,Probe_test] = CalcSensFromDesign(round1.ProbeOpt(2).Coil.N_v,round1.ProbeOpt(2).Coil.N_l,round1.ProbeOpt(2).Coil.do,round1.ProbeOpt(2).Coil.di,1) ; % (NrLoopsPerLayer,NrLayers,OutWireDiam,CoopDiam, AmplMode)