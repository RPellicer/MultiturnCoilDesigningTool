%% Magnetic model Far field
% emf_far=  ?B ?_(i=1:N_t)[?*r_i^2]

function emf_far = emfFar(freqs,r)
    w = 2*pi*freqs;
    Seq = areaTot(r);    % Calculate the total area of the coil
    emf_far = 1j*w.*Seq; % per B unity!
    %     figure,loglog(freqtmp2,abs(emf_far )),grid on, title('Field transfer coefficient')
    %     hold on
    %     plot(3300,70e3,'*r')
    %     legend('Calculated', 'Reported by Matlashov')
end
