clc
%% This script shows litz wire calculations
% Define wire properties
Coil2_5 = struct('diam',2.5,'diamOut',2.7,'area',4.909,'resistance',3.46e-3,'mass',44.424);
Coil0_8 = struct('diam',0.8,'diamOut',0.9,'area',0.503,'resistance',3.38e-2,'mass',4.549);

equivCoilNr_resistance = Coil0_8.resistance / Coil2_5.resistance;
integCoilNr = ceil(equivCoilNr_resistance);
equivCoilNr_area = Coil2_5.area / Coil0_8.area;
equivCoilNr_mass = Coil2_5.mass / Coil0_8.mass;

diffCoilwinding_diam = Coil0_8.diamOut*4/Coil2_5.diamOut; % See Bp-HomeMade10litzWire.pptx
diffCoilNr_resistance = Coil0_8.resistance / (Coil2_5.resistance*integCoilNr);
diffCoilNr_area = (Coil0_8.area*integCoilNr)/Coil2_5.area;
diffCoilNr_mass = (Coil0_8.mass*integCoilNr)/Coil2_5.mass;


disp(['To keep the resistance of current coil we need ',num2str(equivCoilNr_resistance),' wires of diameter ', num2str(Coil0_8.diam)])
disp(['This leads to a mass difference of ',num2str(diffCoilNr_mass),' (extra material)'])
disp(['And a resistance difference of ',num2str(diffCoilNr_mass),' less resistance'])
disp(['Needing ',num2str(diffCoilNr_area),' more space to be winded'])

%% Estimate actual coil's length
R2_5 = 0.2577; % Ohms
L = R2_5 / Coil2_5.resistance;
disp(['Current coil is about ',num2str(L),' meters long'])
