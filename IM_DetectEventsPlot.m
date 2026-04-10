function IM_DetectEventsPlot(Event)
% ---- time conversion
sec2day=1/60/60/24;
% sec2year=sec2day;
sec2year=sec2day/365;

% ==============================================
% fig=figure;
% tiledlayout('flow','Padding','compact','TileSpacing','loose');
% nexttile; hold on;
Fun_defaultAxes;
eTime=Event.Time*sec2year;
eMag=Event.Mw;
stem(eTime,eMag,'k','Marker','o','MarkerEdgeColor','k','LineWidth',0.6,'MarkerSize',14, ...
    'BaseValue',round(min(eMag)-1));
xlim([0.98*min(eTime) 1.02*max(eTime)])
ylim([round(min(eMag)-1) round(max(eMag)+1)]);
ylabel('Magnitude');
xlabel('Time');


end
