% Plot the effect of the Q value of the capacitors
M = [6.9608e-15
5.2168e-15
4.4679e-15
4.3679e-15
4.1123e-15]
% 7.0997e-15
% 5.2733e-15
% 4.4768e-15
% 4.2312e-15

P = [6.6768e-15
4.8126e-15
3.9725e-15
3.8796e-15
3.5256e-15];
% 6.8638e-15
% 4.8609e-15
% 4.0139e-15
% 3.6742e-15];

sweepOptions = [33 100 330 1000 1000000000];


loglog(sweepOptions,M,'-.*r')
hold on, grid on
loglog(sweepOptions,P,'-.*b')
hold off
xlabel('Quality factor of the tunning capacitors')
ylabel('Sensitivity (T/sqrt(Hz))')
legend('Mean', 'Peak')
title(sprintf('Effects of the quality factor of the capacitors \non the sensitivity of the tuned magnetometer'))