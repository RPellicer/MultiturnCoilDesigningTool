function plot_wires(rc,zc,di,do)

circles(squeeze(ones(length(zc),1)*(do/2)), [zc',rc'],'color','b');
hold on
plot(0,0,'*r'), xlabel('Z (m)'), ylabel('r (m)')
% plot(0,0,'*r'), xlabel(['Z (m) nv=' num2str(nv)]), ylabel(['r (m) nl=' num2str(nl)])
circles(squeeze(ones(length(zc),1)*(di/2)), [zc',rc'],'color','g');
plot(zc',rc','+r');
hold on
axis equal, grid on
hold off

end