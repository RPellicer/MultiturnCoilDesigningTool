%% Analyzing efficiency (power losses)(Q = around 100 in the kHz regime) for an inverted L-network
% Han-2006-Analysis and Design of High Efficiency Matching Networks.pdf.
% See also "Insertion loss". pag 16 of lect7_match.pdf (https://www.google.com.au/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&cad=rja&uact=8&ved=0ahUKEwieq-Tao5fNAhXEW5QKHXKtAEgQFggdMAA&url=http%3A%2F%2Frfic.eecs.berkeley.edu%2F142%2Fpdf%2Flect7_match.pdf&usg=AFQjCNF-SM1zdrBem8kFllXPBTCAM1pQkg)

% Rs = Resistance of the source
% Rmatch = Resistance wanted to be seen from the load
% Qc1 = Quality factor of the series capacitor
% Qc2 = Quality factor of the parallel capacitor

function [efficiency, G_effective, G_ideal] = L_inv_effic(Q,Qc1,Qc2)
    efficiency = (1-Q/Qc2)/(1+Q/Qc1); % eqs 16 & 18 from Han-2006-Analysis and Design of High Efficiency Matching Networks.pdf.
    % for Q/Qc << 1 (high efficiency network) it can be simplified to
    % efficiency2 = 1-Q/Qc2-Q/Qc1; % eqs 16 & 18 from Han-2006-Analysis and Design of High Efficiency Matching Networks.pdf. This agrees with lect7_match.pdf
    G_ideal = sqrt(Q^2+1); % from Han-2006-Analysis and Design of High Efficiency Matching Networks.pdf. eq.24
    G_effective = G_ideal*sqrt(efficiency);
end