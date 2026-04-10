function IM_MappingFault2(ITime,IDt,IDisp,IV,IFriction,ITheta,INormalStress,IAccel,IFarDisp,Event,CumEvent)
% =========================================================================
global Fseg Fcen Flen Fangle
global fseg fcen fID fangle flen fcount
global Vini seed
% =========================================================================
% >> Plot parameters:
col=IM_SeedColors(1e3,'bright',seed); % Radom colors for multiplots
% col=Fun_Mycolors2025(1e3);
% col=[0 0 0;1 0 0;0 1 0;0 0 1]; % black-red-green-blue

% windows=500:2500;
windows=1:length(ITime);

Year=60*60*24*365;
KM=1e3;

% =========================================================================
figure;
tiledlayout('flow','TileSpacing','loose','Padding','compact');
nexttile(1);hold on;
% FaultX1=fcen(:,1)+flen(:)/2.*sind(fangle(:));
% FaultX2=fcen(:,1)-flen(:)/2.*sind(fangle(:));
% FaultY1=fcen(:,2)+flen(:)/2.*cosd(fangle(:));
% FaultY2=fcen(:,2)-flen(:)/2.*cosd(fangle(:));
% for i=1:fcount
%     line([FaultX1(i)./1e3,FaultX2(i)./1e3],[FaultY1(i)./1e3,FaultY2(i)./1e3],...
%           'Color','k','linewidth',1.5);
% end
FaultX1=Fcen(:,1)+Flen(:)/2.*sind(Fangle(:));
FaultX2=Fcen(:,1)-Flen(:)/2.*sind(Fangle(:));
FaultY1=Fcen(:,2)+Flen(:)/2.*cosd(Fangle(:));
FaultY2=Fcen(:,2)-Flen(:)/2.*cosd(Fangle(:));
for i=1:Fseg
    line([FaultX1(i)./KM,FaultX2(i)./KM],[FaultY1(i)./KM,FaultY2(i)./KM], ...
        'Color',col(i,:),'linewidth',2.5);
end
% scatter(fcen(:,1)./KM,fcen(:,2)./KM,40,'MarkerFaceColor','k','MarkerEdgeColor', ...
%     'k','Marker','x','LineWidth',1.5); % plot segments
if Event.Count~=0
    Occur=Event.Location;
    Ox=Occur(:,1);
    Oy=Occur(:,2);
    % scatter(Ox./KM,Oy./KM,150,Event.Time/Year,"filled",'MarkerEdgeColor', ...
    %     'k','Marker','p','LineWidth',0.6); % with color coding
    colormap(flip(slanCM('rainbow',20)));
    scatter(Ox./KM,Oy./KM,60,'k',"filled",'MarkerEdgeColor', ...
        'k','Marker','p','LineWidth',0.6);
end
Fun_defaultAxes;grid;
set(gca,'color',[.97 .97 .97]);
xlabel('XDist');
ylabel('YDist');


% nexttile(2);
% plot(fID,Vini,'LineWidth',2,'Color','k','Marker','o','MarkerSize',12,'MarkerFaceColor','k');
% ax=gca;
% ax.YScale="log";
% Fun_defaultAxes;grid;
% xlabel('Fault ID');
% ylabel('Initial Velocity');
nexttile(2);hold on;
if Event.Count~=0
    IM_DetectEventsPlot(ITime,Event,CumEvent);
end

nexttile(3);hold on;
plot(ITime./Year,log10(real(IV())),"LineWidth",1.8);
Fun_defaultAxes;ylabel('Velocity');xlabel('Time');


nexttile(4);hold on;
plot(ITime./Year,log10(real(IDisp())),"LineWidth",1.8);
Fun_defaultAxes;ylabel('Displacement');xlabel('Time');

set(gcf,'position',[200,50,1100,820]);


end