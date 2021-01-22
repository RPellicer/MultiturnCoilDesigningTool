Inoise = 1;
Vnoise = 1;
R1 = logspace(-6,6,100);
Rf = 5000;
G = Rf./R1;
Reff = R1.*Rf./(R1+Rf);
VInOut = G.*(Inoise.*Reff);
VVnOut = G.*Vnoise;
VOutNoiseTot = sqrt(VInOut.^2+VVnOut.^2);
VInin = (Inoise.*Reff);
VVnin = (Vnoise.*G)./G;
VinNoiseTot = sqrt(VInin.^2+VVnin.^2);
loglog(R1,VInOut)
hold on
loglog(R1,VVnOut)
loglog(R1,VOutNoiseTot)
legend('In','Vn','Tot')
hold off
h=figure;
loglog(R1,VInin)
hold on
loglog(R1,VVnin)
loglog(R1,VinNoiseTot)
legend('In','Vn','Tot')
hold off