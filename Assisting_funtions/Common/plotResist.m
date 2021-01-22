function plotResist(R_dc, R_ac, F, G_2, freq)
subplot(2,1,1)
semilogx(freq,R_ac,'-*'), ylabel('resistance(ohm)'),xlabel('frequency(Hz)')
hold on
semilogx(freq,F*R_dc,'r')
semilogx(freq,G_2*R_dc,'g')
semilogx(freq,R_dc,'k')
legend('Total Resistance(AC)','Skin depth EFFECT','Proximity EFFECT','DC resistance')
hold off
subplot(2,1,2)
loglog(freq,R_ac,'-*'), ylabel('resistance(ohm)'),xlabel('frequency(Hz)')
hold on
loglog(freq,F*R_dc,'r')
loglog(freq,G_2*R_dc,'g')
loglog(freq,R_dc,'-*k')
legend('Total Resistance(AC)','Skin depth EFFECT','Proximity EFFECT','DC resistance')
grid on
hold off
end