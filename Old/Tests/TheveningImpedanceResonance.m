clear cll, close all, clc
Rb = 71.7;
Lb = 57.4e-3;
Cb = 8.2e-12;
f = logspace(0,9,1000);
w = 2*pi*f;
Zth = (Rb+1i.*w.*Lb)./((1-Lb.*Cb.*w.^2)+1i.*w.*Rb.*Cb);
figure(1),
subplot(3,1,1),
semilogx(f,real(Zth));
title('real')
subplot(3,1,2),
semilogx(f,imag(Zth));
title('imaginary')
subplot(3,1,3),
semilogx(f,abs(Zth));
title('magnitude')