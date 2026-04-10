function [Nrate,Time2inst]=IM_calAnnualSeismicity(Event,ITime,CumEvent)
% =========================================================================
global a H Vini TMax TMin fseg fID
Fsegeff=length(fseg);
Year=60*60*24*365;
% Yealy seismicity rate:
Nrate=Fsegeff/((TMax-TMin)/Year);

% Find the theoretical time (Time2inst) when rupture should occur on each
% fault and the actual rupture time (Event.Time)
% Given that nearly no tectonic loading, the instability time is:
Time2inst=a./(H*Vini);
% =========================================================================
commonID=intersect(unique(fID(:)), unique(Event.FaultID(:)));
maskTheory=ismember(fID(:), commonID);
maskEvent=ismember(Event.FaultID(:), commonID);
figure;Fun_defaultAxes; hold on;
stairs(fID(maskTheory),Time2inst(maskTheory)/Year,'LineWidth',0.5,'Color','k','Marker','o', ...
    'MarkerSize',8,'MarkerFaceColor',[.3 .3 .3],'DisplayName','Theory');
stairs(Event.FaultID(maskEvent), Event.Time(maskEvent)/Year,'LineWidth',0.6,'Color','r', ...
    'Marker','s','MarkerSize',10,'MarkerFaceColor','none','DisplayName','Events');
legend('Location','southeast','FontSize',18);
ax=gca;
ax.YScale="log";
xlabel('Fault ID');
ylabel('Time to inst');
set(gcf,'position',[300,50,1000,580]);

figure;Fun_defaultAxes; hold on;
plot(ITime/Year,CumEvent,'LineWidth',1.8,'Color','r'); % True
seq=linspace(0,max(Event.Time),10);
plot(seq./Year,Nrate*seq./Year,'LineWidth',2.4,'Color',[0.3373    0.7882    0.0157]); % Annual
% Adding the theory
if min(Time2inst)<max(Event.Time)
    xline(Time2inst(Time2inst<max(Event.Time))/Year,'LineWidth',1,'Color',[.5 .5 .5],'LineStyle',':');
end
set(gcf,'position',[300,50,700,650]);
xlabel('Time');
ylabel('Cum.Events');

end