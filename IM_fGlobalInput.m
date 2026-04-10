function IM_fGlobalInput()
% =========================================================================
% Common globals:
global Doing
global seed G nu rho Width TotalStep TotalRecord Split Dt CCrit VCrit ACrit MomentCrit MagCrit
global a b Dc V0 VlStress k0 Normal0 H CritLen
global Fseg Fcen Flen Fangle FRL FiniNormal FiniShear  % Fault-params
global fcen flen fangle fRL finiNormal finiShear fcount fseg flensum fID  % Segment-params
global TMax TMin Fv
% >> Seed parameters
seed=2025; % seed=randi([1 1e3],1);
rng(seed);
% =========================================================================
% >> Property parameters
G=20e9;    % Shear modulus (pa)
nu=0.2;    % Poisson ratio
rho=2400;  % Rock density (kg/m^3)
Width=800; % Fault width (m)
CritLen=1000;     % Fault Segment length (m)
% CritLen=2*G*Dc/(pi*InitialSigma(1)*(b-a));
MomentCrit=1e10; % Occurrence Threshold (Nm)
MagCrit=2;       % Detection Threshold

% =========================================================================
% >> Process parameters
TotalStep=1.5e5;
Split=20;
TotalRecord=TotalStep/Split;
Dt=1;       % Self-adjusted time-step (any >1e-5 s)

% =========================================================================
% >> Friction parameters
CCrit=1e-8; % Convergence Critirion
VCrit=1e-6; % Threshold that may cause acceleration behavior (m/s)
ACrit=1e-3; % Threshold similar to VCrit (m/s^2)
a=0.003;    % RSF params.a
b=0.006;    % RSF params.b
Dc=100e-6;  % 1e-6 (m)~1 (μm)
V0=1e-9;   % Reference velocity
k0=1e8;             % Fault stiffness (used in H)
Normal0=15e6;       % Fault normal stress (used in H)
H=-k0/Normal0+b/Dc; % H constant

% =========================================================================
% >> Velocity parameters
% ======= Method <1> =======
% From the formulation when tectonic loading is week:
% The time to instability: t=A/(H*Vi)
% And we can obtain the expected yealy seismicity: N(faults)/(TMax-TMin)
TMax=2e9;  % (s)
TMin=1e1;  % (s)
Time2inst=TMin+rand(Fseg,1)*(TMax-TMin); % Random
% TRand=TMin+0.2*ones(Fseg,1)*(TMax-TMin); % Same
% TRand=sort(TMin+rand(Fseg,1)*(TMax-TMin)); % Main-Afters
Fv=a./(H.*Time2inst);

% ======= Method <2> =======
% Fv=rand(Fseg,1).^3*1e-13;   % Radom
% Fv=sort(rand(Fseg,1).^3*1e-13,'descend');   % Main-Afters


% =========================================================================
% >> Fault parameters 1
if Doing==1
    Fseg=4;
    Fcen=[
        6000,6000;
        2000,5000;
        10000,2000;
        8000,10000];
    Flen=[5;1;1;1]*1e3;
    Fangle=[30;30;30;30];   % positive to the right
    FRL=[-1;-1;-1;-1];      % positive: right-lateral
    VlStress=ones(Fseg,1)*cosd(30)*1e-6;     % Far-field-Loading Rate (Pa/s;e.g. 1e-3)
    % Notice value of cosd(30) is related to Fangle and will extend dimension
    % next step
    % VlStress=[0 0 1e-2 0]';            % Far-field-Loading Rate
    FiniNormal=ones(Fseg,1)*15e6;   % Initial normal stress
    FiniShear=FiniNormal*0.6;       % Initial shear stress
    [fcen,flen,fangle,~,fRL,finiNormal,finiShear,VlStress,fcount,fseg,flensum,fID]...
        =IM_DiscrestFault();        % Only Discrete
end
% =========================================================================
% >> Fault parameters 2
if Doing==2 % Needs more parameters
    global Fangle1 Fangle2 FlenMin FlenMax DxMin DxMax DyMin DyMax
    global Sigma1 Sigma3 MainPart
    Fseg=500;
    Fangle1=30;
    Fangle2=30;
    FlenMin=100;
    FlenMax=5*1e3;
    DxMin=0;
    DxMax=40*1e3;
    DyMin=0;
    DyMax=40*1e3;
    Sigma1=30e6;    % Maximum principal stress (also related to Fangle)
    Sigma3=10e6;    % Minimum principal stress (also related to Fangle)
    MainPart=0.3;   % Ratio of Main-Fault (%)
    [fcen, flen, fangle, ~, fRL, finiNormal, finiShear, fcount, ...
    fseg, flensum, fID, Fcen, Flen, Fangle, ~, FRL, FiniNormal, FiniShear]...
        =IM_FaultGenerate();  % Generate and Discrete
    VlStress=ones(fcount,1)*cosd(30)*1e-6;    % Far-field-Loading Rate (Pa/s)
end
% =========================================================================
fprintf('# Global Data have been loaded\n');

