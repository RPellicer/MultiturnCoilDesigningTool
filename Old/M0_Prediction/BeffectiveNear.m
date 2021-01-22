%% The amplitude of the magnetic field produced by the magnetization of the voxel
% from "Calculated signal-to-noise ratio of MRI detected with SQUIDs and
% Faraday detectors in fields from 10 lT to 1.5 T"

function Bdet = BeffectiveNear(M,Vvoxel,Atot,BrecipPerp)
    u0 = 4*pi*1e-7; %permeability of free space
    Bdet = u0*BrecipPerp*M*Vvoxel/(4*pi*Atot);
end