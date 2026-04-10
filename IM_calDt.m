function [DtRef]=IM_calDt(V,SlowOrFast)
% Based on the current maximum values of velocity, provide a reasonable time-step
% empirically
% =========================================================================
% If the system is in a rapid or potentially unstable phase, adjust the 
% time step to a smaller value
if SlowOrFast==1
    DtLower=5e-2; % Time step [second]
    DtUpper=1e-2; % Time step [second]
    VLower=1e-4;
    VUpper=1e-2;
    Vmax=max(V);
    
    if Vmax<VLower
        DtRef=DtLower;
    elseif Vmax > VLower && Vmax < VUpper
        x=(log10(Vmax)-log10(VLower))/(log10(VUpper)-log10(VLower));
        y=3*x^2-2*x^3;
        LogDt=log10(DtLower)+y*(log10(DtUpper)-log10(DtLower));
        DtRef=10^LogDt;
        %     DtRef=Dt_Lower-(Dt_Lower-Dt_Upper)*(3*((Vmax-VTreshLower)/(VTreshUpper-VTreshLower))^2-2*((Vmax-VTreshLower)/(VTreshUpper-VTreshLower))^3);
    else
        DtRef=DtUpper;
    end
% =========================================================================
else
    % The larger the value of <DtLower>, the greater the total duration or 
    % the longer the creep time
    % 3600 s correspond to 1 h
    DtLower=3600; % Time step [second]
    DtUpper=5e-2; % Time step [second]
    VLower=1e-7;
    VUpper=1e-5;
    Vmax=max(V);

    if Vmax<VLower
        DtRef=DtLower;
    elseif Vmax > VLower && Vmax < VUpper
        x=(log10(Vmax)-log10(VLower))/(log10(VUpper)-log10(VLower));
        y=3*x^2-2*x^3;
        LogDt=log10(DtLower)+y*(log10(DtUpper)-log10(DtLower));
        DtRef=10^LogDt;
        % DtRef=Dt_Lower-(Dt_Lower-Dt_Upper)*(3*((Vmax-VTreshLower)/
        %     (VTreshUpper-VTreshLower))^2-2*((Vmax-VTreshLower)/(VTreshUpper-VTreshLower))^3);
    else
        DtRef=DtUpper;
    end    
    
end


end
