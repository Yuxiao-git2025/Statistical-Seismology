function IM_MappingFault1(ITime,IDt,IDisp,IV,IFriction,ITheta,INormalStress,IAccel,IFarDisp,Event)
% =========================================================================
global Fcen Flen Fangle
global fcen flen fangle
global seed MagCrit FnumMain fnumMain FnumBkgd fnumBkgd
% global fseg Fseg
% =========================================================================
% >> Plot parameters:
rng(seed);
% col=IM_SeedColors(1e3,'bright',seed); % Radom colors for multiplots
% col=Fun_Mycolors2025(1e3);
% col=[0 0 0;1 0 0;0 1 0;0 0 1]; % black-red-green-blue for example

Year=60*60*24*365;
KM=1e3;
% Notice: True faults number could be derivated from length(fseg) !
Fsegeff=FnumMain+FnumBkgd;
fsegeff=fnumMain+fnumBkgd;
fprintf('# After removal\n');
fprintf('# The total %d faults is divided into: %d (Main) + %d (Secondary)\n',Fsegeff,FnumMain,FnumBkgd);
fprintf('# The total %d segments is divided into: %d (Main) + %d (Secondary)\n',fsegeff,fnumMain,fnumBkgd);
% =========================================================================

figure;
tiledlayout(1,2,'TileSpacing','loose','Padding','compact');
% Faults Mapping:
nexttile(1);hold on;
% FaultX1=Fcen(:,1)+Flen(:)/2.*sind(Fangle(:));
% FaultX2=Fcen(:,1)-Flen(:)/2.*sind(Fangle(:));
% FaultY1=Fcen(:,2)+Flen(:)/2.*cosd(Fangle(:));
% FaultY2=Fcen(:,2)-Flen(:)/2.*cosd(Fangle(:));
% Just plot faults after removed
FaultX1=fcen(:,1)+flen(:)/2.*sind(fangle(:));
FaultX2=fcen(:,1)-flen(:)/2.*sind(fangle(:));
FaultY1=fcen(:,2)+flen(:)/2.*cosd(fangle(:));
FaultY2=fcen(:,2)-flen(:)/2.*cosd(fangle(:));
for i=1:fnumMain
    line([FaultX1(i)./KM,FaultX2(i)./KM],[FaultY1(i)./KM,FaultY2(i)./KM], ...
        'Color','r','linewidth',1.6);
end
for i=fnumMain+1:fnumMain+fnumBkgd
    line([FaultX1(i)./KM,FaultX2(i)./KM],[FaultY1(i)./KM,FaultY2(i)./KM], ...
        'Color','k','linewidth',1.2);
end
% Event Mapping:
if Event.Count~=0
    Eloc=Event.Location;
    Ex=Eloc(:,1);
    Ey=Eloc(:,2);
    scatter(Ex./KM,Ey./KM,8.^(Event.Mw-MagCrit+0.1)+20,Event.Time/Year,"filled",'MarkerEdgeColor', ...
        'flat','Marker','o','LineWidth',0.6); % with color coding
    colormap((slanCM('rainbow',12)));
    c=colorbar;
    c.LineWidth=1; c.Label.String='Time';
%     scatter(Ox./KM,Oy./KM,150,'k','MarkerEdgeColor', ...
%         'k','Marker','p','LineWidth',0.6);
end
Fun_defaultAxes;grid;
set(gca,'color',[.97 .97 .97]);
xlabel('XDist');
ylabel('YDist');

% M-T plot:
nexttile(2);hold on;
if Event.Count~=0
    IM_DetectEventsPlot(Event);
end
set(gcf,'position',[200,120,1200,520]);


% =========================================================================
% >> Other plots

% nexttile(3);hold on;
% p1=1;
% p2=fseg(1);
% for i=1:Fsegeff
%     pid=p1:p2;
%     if i~=Fsegeff
%         p1=p2+1;
%         p2=p2+fseg(i+1);
%     end
%     plot(ITime(:)./Year,log10(real(IV(:,pid(1)))), ...
%         "LineWidth",1.8,'Color',col(i,:));
% end
% Fun_defaultAxes;ylabel('Velocity');xlabel('Time');

% figure;
% global fID Vini
% plot(fID,Vini,'LineWidth',1.2,'Color','k','Marker','o','MarkerSize',6,'MarkerFaceColor','k');
% ax=gca;
% ax.YScale="log";
% Fun_defaultAxes;grid;
% xlabel('Fault ID');
% ylabel('Initial Velocity');

% =========================================================================
% Location mapping:
% figure;
% tiledlayout(1,2,'TileSpacing','loose','Padding','compact');
% nexttile;
% scatter(Event.Time./Year,Event.Location(:,2)./KM,80,'MarkerFaceColor',[0.8706    0.8706    0.8706], ...
%     'MarkerEdgeColor','k');
% Fun_defaultAxes;ylabel('YDist');xlabel('Time');
% 
% nexttile;
% scatter(Event.Time./Year,Event.Location(:,1)./KM,80,'MarkerFaceColor',[0.8706    0.8706    0.8706], ...
%     'MarkerEdgeColor','k');
% Fun_defaultAxes;ylabel('XDist');xlabel('Time');
% set(gcf,'position',[200,50,1200,520]);
end