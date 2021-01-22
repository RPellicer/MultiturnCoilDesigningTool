% Calculate an estimated and normalised frequency gain curve by the 
% freqRespLinv(71.7,57.4e-3,4.8629e-9,4.5584e-8,2e-12,60e6,[100 1000])
function gain = freqRespLinv(Rs,L,C1,C2,Ca,Ra,freqs)
    s = 1j*2*pi*freqs;
    gain = -(C1*Ra.*s)./(Ca*Ra.*s + C1*Ra.*s + C2*Ra.*s + C1*Rs.*s + C1*L.*s.^2 +...
    Ca*C1*L*Ra.*s.^3 + C1*C2*L*Ra.*s.^3 + Ca*C1*Ra*Rs.*s.^2 + C1*C2*Ra*Rs.*s.^2 + 1);

%     figure(20)
%     plot(freqs,abs(gain)), grid on
end