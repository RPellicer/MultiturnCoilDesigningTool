% Create an estimated and normalised frequency gain curve by the 
% Use only when gain function is not known because it can change
% significantly the results
% #################################
%          witch of Agnesi
% #################################
% is of importance to physicists because it approximates the spectral energy...
% distribution of x-ray lines and optical lines, as well as the power ...
% dissipated in sharply tuned resonant circuits. The maximum height is h...
% while the half-width at half-maximum is a. Henceforth, we shall designate a as the “half-width.”
%         y=h*(a^2)./(a^2+x.^2);
% freqResp(1e4,20,9000:0.1:11000)

function gain = freqRespWitchAgnesi(f0,Q,freqs)
    x = (freqs-f0)/f0;
    h = 1;
    a = 1/2;
    gain=(h*(a^2)./(a^2+(Q*x).^2));
%     figure(10)
%     plot(freqs,gain), grid on
end