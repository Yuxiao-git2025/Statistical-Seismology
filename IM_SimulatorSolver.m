function [ITime,IDt,IDisp,IV,IFriction,ITheta,INormalStress,IAccel,IFarDisp]=IM_SimulatorSolver()
% =========================================================================
global seed G nu rho TotalStep TotalRecord Split Dt CCrit VCrit ACrit 
global a b Dc V0 VlStress
global Vini Thetaini kShear kNormal EShear ENormal
global Fv
global fcen flen fangle fRL finiNormal finiShear fcount fseg flensum fID
rng(seed);

% =========================================================================
% % >> Discretize the faults (Now have finished in Parameters-input function)
% [fcen,flen,fangle,frad,fRL,finiNormal,finiShear,VlStress,fcount,fseg,flensum,fID]...
%     =IM_DiscrestFault();

% >> Generate Stiffness-Matrix (kij) and k1, k2 for Shear and Normal Matrix
[kShear,kNormal]=IM_calStiffness(fcount,fcen,flen,fangle,G,nu,fRL);
% Remove faults that too-close
count1=fcount;
[fcen,flen,fangle,fRL,finiNormal,finiShear,fcount,flensum,fID,VlStress,fseg]=...
    IM_calCloseness(kShear,kNormal);
count2=fcount;
% Recalculate the Stiffness
[kShear,kNormal]=IM_calStiffness(fcount,fcen,flen,fangle,G,nu,fRL);
clc;
datime1=datetime('now');
fprintf('# Have Removed %d segments\n',count1-count2);
fprintf('# Total steps is %.0e and every %d steps data will be saved\n',TotalStep,Split);
fprintf('# Have Created %d number of faults and %d segments\n',length(fseg),sum(fseg));

% =========================================================================
% >> Initialize velocity for each segment
% Notice: 
% When the velocity magnitude is relatively large (i.e. >1e-12),
% instability will occur rapidly, and it is possible that secondary faults
% are more prone to rupture than the main fault (because they are weaker in
% strength); When the speed is relatively low, the time to instability is
% longer
Vini=zeros(fcount,1);
if max(fID)>length(Fv); error('# Please check velocity input'); end
for i=1:fcount
    Vini(i)=Fv(fID(i));
end

% >> Expand the parameters
% The larger the status value is, the better the initial healing is,
% The time of instability will be advanced
Thetaini=ones(fcount,1)*1e10; 
a=a*ones(fcount,1);
b=b*ones(fcount,1);
Dc=Dc*ones(fcount,1);

% >> Calculate Normalized Shear- and Normal- stiffness (dimensionless)
% Represent the relative coupling shear/normal strength of each fault segment
% It won't be used ENormal here
EShear=zeros(fcount,fcount);
ENormal=zeros(fcount,fcount);
Omega=zeros(fcount,1);
kShearii=zeros(fcount,1);
alpha=1;
Mass=(alpha*rho*flensum)/((1-nu)*pi^2); % Mass per unit contact area
% Notice: L is flensum but not flen
for i=1:fcount
    kShearii(i)=kShear(i,i);
    Omega(i)=sqrt(kShearii(i)/Mass(i));
    for j=1:fcount
        if i==j
            EShear(i,j)=0;
        else
            EShear(i,j)=kShear(i,j)/kShear(i,i);
            ENormal(i,j)=kNormal(i,j)/kShear(i,i);
        end
    end
end
FrictionI=finiShear./finiNormal; % initial friction each segment
Friction0=FrictionI-a.*log(Vini./V0)-b.*log(Thetaini.*V0./Dc); % Reference friction
for i=1:length(Vini)
    if (Friction0(i)<=0) || (Friction0(i)>1)
%         fprintf('# Check the value of mu0 \n');
        Friction0(i)=0.6; % for those v=0
    end
end

% >> Data initialize
ITime=zeros(TotalRecord,1);
IDisp=zeros(TotalRecord,fcount);
IV=zeros(TotalRecord,fcount);
IFriction=zeros(TotalRecord,fcount);
ITheta=zeros(TotalRecord,fcount);
INormalStress=zeros(TotalRecord,fcount);
IAccel=zeros(TotalRecord,fcount);
IDt=zeros(TotalRecord,1);
IFarDisp=zeros(TotalRecord,fcount);
% >> Step initialize
% Given that: Ks*u=τ, i.e. u=τ/Ks=Ks\τ
iniDisp=(kShear\finiShear); % initial load point
FarDisp_Old=iniDisp;
Disp=zeros(fcount,1);
DispOld=zeros(fcount,1);

Accel=0;
TOld=0;
T=0;

Friction=FrictionI;
FrictionOld=Friction;

Instability=zeros(fcount,1);

AccelOld=zeros(fcount,1);

Theta=Thetaini;
ThetaOld=Theta;

VOld=Vini;
Vlp=VlStress./kShearii;

EffSigma_Old=finiNormal;
EffSigma=zeros(fcount,1);
IsInstability=zeros(fcount,1);

SlowOrFast=0; % initial is slow velocity that corresponds to a larger step
Solver=zeros(fcount,1); % 0 for slow and 1 for rapid
Step=0;
% =========================================================================
% >> Loop beings
% =========================================================================
hbar=waitbar(0,'Processing','Name','Processing Starting');
for i=1:TotalStep
    V=VOld;
    [DtRef]=IM_calDt(V,SlowOrFast);
    if Dt<DtRef/100
        Dt=Dt*100; % Gradully increase to DtRef
    else
        Dt=DtRef;  % Given a limit value
    end
    if max(Instability)==2
        Dt=Dt*100;
    end

    Break=0;
    while Break==0
        DtOld=Dt;
        Disp_i=DispOld;
        Friction_i=FrictionOld;
        V_i=VOld;
        Theta_i=zeros(fcount,1);
        % Assign Dt to each segment
        Dtseg=ones(fcount,1)*Dt;
        % velocity perturbation
        DV=10.^(log10(VOld)-(log10(VOld-AccelOld*DtOld)-log10(VOld))*Dt/DtOld);

        EffSigma=finiNormal+kNormal*DispOld;
        for k=1:fcount
            if EffSigma(k)<1e6
                EffSigma(k)=1e6; % Minimum
            end
        end
        % Displacement of the far-field loading point varies with time:
        FarDisp=FarDisp_Old+Vlp.*Dt;
        % Considering the elastic interaction between segments:
        EDisp=EShear*(FarDisp-DispOld); % dδ=δlp-δ
        EVel=EShear*(Vlp-VOld);
        % Total disp:
        TotalDisp=FarDisp+EDisp; 
        
        % parfor FaultIdx=1:FaultCount
        for id=1:fcount
            if Solver(id)==1 % case of instability
                 [V(id),Friction(id),Disp(id),Theta(id),EffSigma(id),Dtseg(id),IsInstability(id)]...
                    =IM_SmallStep(CCrit,DispOld(id),FrictionOld(id),ThetaOld(id),VOld(id),...
                    Friction0(id),a(id),b(id),Dc(id),V0,kShearii(id),Dt,Omega(id),G, Instability(id),...
                    Disp_i(id), Friction_i(id), V_i(id), Theta_i(id), DV(id), EffSigma(id),TotalDisp(id));
            else             % case of creep/loading
                [V(id),Friction(id),Disp(id),Theta(id),EffSigma(id),Dtseg(id),IsInstability(id)]...
                    =IM_LargeStep(CCrit,DispOld(id),FrictionOld(id),ThetaOld(id),VOld(id),...
                    Friction0(id),a(id),b(id),Dc(id),V0,kShearii(id),Dt,Vlp(id),AccelOld(id),...
                    Mass(id),G, EffSigma_Old(id),Instability(id),Disp_i(id), ...
                    Friction_i(id),V_i(id),Theta_i(id),DV(id),EVel(id),EffSigma(id));
            end
            if IsInstability(id)>0
                Solver(id)=1;
            end
        end

        if min(Dtseg)<max(Dtseg)
            if max(IsInstability)==3 % for 0/1/3, 3 is most serious
                Dt=1e-5;
            else
                Dt=min(Dtseg);
            end
        else
            Break=1; % Break if all segments are accepted same Dt
        end
    end
    Dt=min(Dtseg);

    % >> Data save
    if rem(i,Split)==0       
        Step=Step+1;
        ITime(Step,1)=T;
        IDt(Step,1)=Dt;
        IDisp(Step,:)=Disp;
        IV(Step,:)=V;
        IFriction(Step,:)=Friction;
        ITheta(Step,:)=Theta;
        INormalStress(Step,:)=EffSigma;
        IAccel(Step,:)=Accel;
        IFarDisp(Step,:)=FarDisp;
    end
    % >> Step save
    T=TOld+Dt; % step T
    TOld=T;
    Accel=(V-VOld)/Dt;
    AccelOld=Accel;
    DispOld=Disp;
    ThetaOld=Theta;
    VOld=V;
    FrictionOld=Friction;
    FarDisp_Old=FarDisp;
    EffSigma_Old=EffSigma;
    Instability=IsInstability;
    
    % >> Determine whether the time step should be adjusted
    for idx=1:fcount
        if V(idx)>VCrit || abs(Accel(idx))>ACrit
            Solver(idx)=1;
        else
            Solver(idx)=0;
        end
    end
    if any(Solver>0)
        SlowOrFast=1;
    else
        SlowOrFast=0;
    end
    if mod(i,TotalStep/20)==0
        waitbar(i/TotalStep, hbar, ['Now ' num2str(i/TotalStep*100) ' %']);
    end
end
% =========================================================================
datime2=datetime('now');
fprintf('# Now Finished  %s\n',datime2);
fprintf('# Using %s times\n',(datime2-datime1));
close(hbar);
% =========================================================================

end