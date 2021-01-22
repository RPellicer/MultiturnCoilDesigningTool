clear all, close all, clc
fname='ModelTunedCircuitTest1.cir';
% fname='ModelSymTunedCircuitTest1.cir';
scam
resp = v_3/Vin;
% convert the expression to a MATLAB transfer function object
% First we separate out the numerator and denominator
% [n,d]=numden(eval(resp));
w = sym('w')
% resp = subs(resp,s,w*1i);
[nsym,dsym]=numden(resp);
[n,d]=numden(eval(resp));
% convert them to MATLAB polynomials
mySys=tf(sym2poly(n),sym2poly(d));
% Step Response
figure(1)
step(mySys)
% Bode Plot
figure(2)
opts = bodeoptions;
opts.FreqUnits = 'Hz';
h = bodeplot(mySys,opts);
grid on
% Find roots
roots = solve(dsym==0,s);
w0 = imag(eval(roots));
f0 = w0./(2*pi);
% Calculate the values at resonance
resp_w0 = eval(subs(resp,s,(w0.*1i)));
% Derivate
dresp = diff(resp,s);
% dresp = diff(resp,w);
[dnsym,ddsym]=numden(dresp);
droots = solve(dnsym==0,s);
dw0 = eval(droots);
df0 = dw0/(2*pi);
dresp_dw0 = eval(subs(resp,s,(dw0.*1i)));

% Tests
omega = 1i*logspace(4,6,10000);
AC_response = eval(subs(resp,s,omega));
figure(3);
semilogx(abs(omega/(2*pi)),10*log10(abs(AC_response)));
grid on