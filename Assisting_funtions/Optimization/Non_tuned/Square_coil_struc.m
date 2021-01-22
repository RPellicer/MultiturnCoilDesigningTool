%% This function assumes the design to be square--> N wires axially = N wires radially
% percentWire --> The percent diameter of the copper wire in total Diameter(Copper + Insulation) 
function TotalSensitivity = Square_coil_struc(Nwires,percentWire) 
    %Define constants
	% the structure of Coil 
    rout = 0.03;
    rin = 0.01;
    freq = 3.3e3;
    %% % Design by number of turns
    do = ((rout -rin)/Nwires); % [METERS] wire outer diameter insulation included
    di = do*(percentWire)/100; % [METERS] cooper outer diameter
    % To match indutance
    N_v = Nwires; % Number of turns per layer
    N_l = Nwires; % Number of layers
    dr = 0.0; % axial filling factor = a-sin(pi/3)*a in an hexagonal setup
    dz = 0.0; % azimut filling factor
%     s = do + sqrt(dr^2+dz^2); % distance between wire centres
    R_i = 0.01; % Coil internal radius

    % Noise from the amplifier (fixed to the amplifier type INA217)
    i_n = 0.8e-12;          %INA217
    e_n = 1.2e-9;           %INA217
    e_OutStagNoise = 90e-9; %INA217
% %     Ideal amplifier
%     i_n = 0.8e-20; %INA217
%     e_n = 1.2e-20;  %INA217
%     e_OutStagNoise = 90e-20;

    Rf = 2*5000; % Feedback resistance of the amplifier (2x in the case if the INA)

    % Frequencies of interest
    freqs= freq;
    % freqs= 3.3e3;

    TotalSensitivity = sensitivityI2V(i_n,e_n,e_OutStagNoise,Rf,N_v,N_l,R_i,di,do,dr,dz,freqs);
   
end
