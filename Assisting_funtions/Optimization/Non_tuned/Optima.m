%% --------------------------------------------------------------------
% This Function assumes the Shape to be Square - layers = coils
Nrlx = 30;  % Number of loops 
% Initialize
sf = zeros(1,Nrlx);
rd = zeros(1,Nrlx);
% Preapare parallel computing
if exist('poolobj')
    delete(poolobj)
end
sfc = [];
poolobj = parpool;
% + 40--> loops are from 1 to 50 --> Net percent is 40% to 90%
for Nloops=1:Nrlx
    s = zeros(1,50);
    parfor RadPerc=1:50
        s(1,RadPerc) = Square_coil_struc(Nloops,(RadPerc+40));
        sfc(Nloops,RadPerc) = s(1,RadPerc);
    end
    [sf(1,Nloops),rd(1,Nloops)] = min(s);
end
% kill the parallel computing
delete(poolobj)

TempLoops=1:Nrlx;
TempRad =(1:50);

figure(2202)
plot (TempLoops,rd+40,'*') 
xlabel('number of wires (square root)') , ylabel('percentage of wire')
grid on

figure(2302)
semilogy(TempLoops,sf)
xlabel('number of wires (square root)') , ylabel('sensitivity')
grid on

figure(2402)
surf(TempRad+40,TempLoops,sfc)
ylabel('number of wires (square root)') , xlabel('percentage of wire'), zlabel('sensitivity')