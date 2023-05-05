function plot_ini_1D(x,xd,Qx,cs0,cvl0,cvr0,ss)
figure(1)
set(gcf,'Renderer', 'painters', 'Position', [100 200 500 500])

subplot(3,1,1)
plot(x,Qx,'b-')
ylabel("$u$ (m s\textsuperscript{-1})","Interpreter","latex")
set(gca,'fontname','times','Layer','top'); box on; grid on;

subplot(3,1,2)
hold on
plot(x,cs0,'k-')
plot(-xd,cvl0,'bs','MarkerFaceColor','b')
plot(+xd,cvr0,'bs','MarkerFaceColor','b')
hold off
ylabel("$c_0$ (unit m\textsuperscript{-3})","Interpreter","latex")
set(gca,'fontname','times','Layer','top'); box on; grid on;

subplot(3,1,3)
plot(x,ss,'k-')
if ss>0; ylim([0 max([ss])*1.1]); end
xlabel("$x$ (m)","Interpreter","latex")
ylabel("$s_s$ (unit m\textsuperscript{-1} s\textsuperscript{-1})","Interpreter","latex")
set(gca,'fontname','times','Layer','top'); box on; grid on;

set(gcf,'color','w'); 
set(gca,'Layer','top')

end