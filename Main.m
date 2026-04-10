% =========================================================================
% >> Create Discrete-Fault-Networks (DFN) and some analysis
% Written by Kyungjae (KJ) Im
% Modified at 2026/04/02-2026/04/08 (Yu)
% =========================================================================
% Many global variables not only include the ones defined by myself, but
% also some that are generated and updated simultaneously during the
% program's execution. Therefore, please analyze the global variables I
% have written with caution
% =========================================================================
% Doing=1 means simple faults-network while Doing=2 means complex-network
global Doing
Doing=1;
IM_fGlobalInput();

%% Basic Calculation
[ITime,IDt,IDisp,IV,IFriction,ITheta,INormalStress,IAccel,IFarDisp]=IM_SimulatorSolver();

%% Detection and Plot
[Event,Fault,CumEvent]=IM_DetectEvents(ITime,IDisp);
IM_MappingFault1(ITime,IDt,IDisp,IV,IFriction,ITheta,INormalStress,IAccel,IFarDisp,Event);

%% More information
% 1. Comparision time-to-instability between theory and true
% 2. Calculate the average seismicity according to faults element
[Nrate,Time2inst]=IM_calAnnualSeismicity(Event,ITime,CumEvent);
