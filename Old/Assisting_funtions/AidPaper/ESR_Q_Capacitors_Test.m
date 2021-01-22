% From http://ksim.kemet.com/Ceramic/CeramicCapSelection.aspx
% Series: C1210C100K8GAC

capESR = [2201000000	113551000	26217000	1615000	512561	57115	21961	11606	5733	2602;
            23139000	1317000	370962	58619	28488	4680	1968	1112	550.505	249.816;
            295075	31203	14582	4833	2621	457.652	194.565	110.734	54.829	24.882;
            9313	2116	1233	473.029	259.852	45.668	19.438	11.075	5.485	2.49;
            729.56	201.565	121.098	47.222	25.978	4.573	1.948	1.113	0.552737	0.251869;
            71.152	20.133	12.138	4.744	2.613	0.46419	0.1995	0.117253	0.05955	0.028066;
            7.299	2.086	1.264	0.497214	0.276798	0.053323	0.024637	0.017647	0.010232	0.005685;
            0.936091	0.28277	0.176617	0.072554	0.043172	0.01224	0.007152	0.007687	0.0053	0.003402];
capFreq = [1,10,100,1000,10000,100000,1000000,10000000];
capVal = [1.00E-11	4.7E-11	1.00E-10	4.7E-10	0.000000001	4.7E-09	0.00000001	0.000000047	0.0000001	0.00000022];
h1 = figure();
surf(capFreq,capVal,capESR',log(capESR'));
set(gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'log')
xlabel('Freq (Hz)'),ylabel('Capacitance (f)'),zlabel('ESR (ohm)')
title('ESR dependence on frequency and capacitance')
cb1 = colorbar;
temp1 = get(cb1,'Ticks');
set(cb1,'TickLabels',round(exp(temp1),2,'significant'));
tempAx1 = h1.CurrentAxes;
axis([tempAx1.XLim tempAx1.YLim min(min(capESR)) max(max(capESR)) min(min(log(capESR))) max(max(log(capESR)))])
view(64,20)

capQ = 1./(2*pi*capFreq'*capVal.*capESR);
h2 = figure();
surf(capFreq,capVal,capQ',log(capQ'));
set(gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'log')
xlabel('Freq (Hz)'),ylabel('Capacitance (f)'),zlabel('Q')
title('Q dependence on frequency and capacitance')
cb2 = colorbar;
temp2 = get(cb2,'Ticks');
set(cb2,'TickLabels',round(exp(temp2),2,'significant'));
tempAx2 = h2.CurrentAxes;
axis([tempAx2.XLim tempAx2.YLim min(min(capQ)) max(max(capQ)) min(min(log(capQ))) max(max(log(capQ)))])
view(64,20)

Reviewer2Prop = 2*pi*capFreq'*capVal.*capESR;
h = figure();
surf(capFreq,capVal,Reviewer2Prop',log(Reviewer2Prop'));
set(gca, 'XScale', 'log', 'YScale', 'log', 'ZScale', 'log')
xlabel('Freq (Hz)'),ylabel('Capacitance (f)'),zlabel('ESR*2*pi*f*C (u)')
title('ESR*2*pi*f*C')