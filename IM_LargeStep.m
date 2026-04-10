function [Return_V,Return_Friction,Return_Disp,Return_Theta,Return_ESigma,Return_Dt,Return_InstabilityThistime]...
    = IM_LargeStep(CCrit,DispOld,FrictionOld,ThetaOld,VOld,...
    Friction0,a,b,Dc,V0,kShearii,Dt,Vl,AccelOld,Mass,G,ESigma_Old,Instability,...
    Disp, Friction, V, Theta, DV, EVel, ESigma)
% Calculate the friction updated by mechanical loading and the friction
% defined by the RSF law
Iteration=0;
IsInstability=0; % 0/1/3
CC=0.1;
CCrit0=CCrit;
% =========================================================================
while CC>CCrit0
    Iteration=Iteration+1;
    VTest=V;

    % >> Original for NR
    Accel=(V-VOld)/Dt;
    % Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law (explicit)
    Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; % Dieterich Law (explicit)

    Friction=FrictionOld+(kShearii/ESigma*(Vl+EVel-V)...
        -(ESigma-ESigma_Old)/Dt*FrictionOld/ESigma...
        -G/(2*3e3)/ESigma*Accel-Mass/ESigma*(Accel-AccelOld)/Dt)*Dt;

    % >> Function Formulation
    % Ultimately, FV is the residual of a nonlinear equation about V
    FV1=V-V0*exp((Friction-Friction0-b*log(V0*Theta/Dc))/a);


    % >> Deviated value for NR
    V=V+DV;
    Accel=(V-VOld)/Dt;
    % Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; Ruina Law
    Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; % Dieterich Law
    Friction=FrictionOld+(kShearii/ESigma*(Vl+EVel-V)...
        -(ESigma-ESigma_Old)/Dt*FrictionOld/ESigma...
        -G/6000/ESigma*Accel-Mass/ESigma*(Accel-AccelOld)/Dt)*Dt;
    FV2=V-V0*exp((Friction-Friction0-b*log(V0*Theta/Dc))/a);

    % >> Update valocity
    V=VTest-FV1/((FV2-FV1)/DV);
    Disp1=DispOld+V*Dt;
    Accel1=(V-VOld)/Dt;
    % =====================================================================
    if V<0
        if Instability==1
            V=1e-20;
            Disp=DispOld+(VOld+V)/2*Dt; % Update disp
            Theta=ThetaOld-V.*ThetaOld./Dc.*log(V.*ThetaOld./Dc)*Dt;
            Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
            IsInstability=1;
            CC=0;
        else
            if  Dt>1e-7
                Dt=Dt/2;
                V=VOld;
                Theta=ThetaOld;
                Friction=FrictionOld;
                Disp=DispOld;
                CC=1;
            else
                V=VOld;%10.^(log10(VOld)-(log10(VOld-Accel*DtOld)-log10(VOld))*Dt/DtOld);
                Instability=1;
                IsInstability=3;
                Disp=DispOld+(VOld+V)/2*Dt; % Update disp
                Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
                Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
                CC=0;
            end
        end
    elseif ~isreal(V)
        if Instability==1
            V=1e-20;
            Disp=DispOld+(VOld+V)/2*Dt; % Update disp
            % Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; % Ruina Law
            Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; % Dieterich Law
            Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
            IsInstability=1;
            CC=0;
        else
            if  Dt>1e-7
                Dt=Dt/2;
                V=VOld;
                Theta=ThetaOld;
                Friction=FrictionOld;
                Disp=DispOld;
                CC=0;
            else
                V=VOld;
                Disp=DispOld+(VOld+V)/2*Dt; % Update disp
                % Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; % Ruina Law
                Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; % Dieterich Law
                Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
                Instability=1;
                IsInstability=1;
                CC=0;
            end
        end

    elseif isnan(V)
        if Instability==1
            V=1e-20;
            Disp=DispOld+(VOld+V)/2*Dt; % Update disp
            % Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; % Ruina Law
            Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; % Dieterich Law
            Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
            IsInstability=1;
            CC=0;
        else
            if  Dt>1e-7
                Dt=Dt/2;
                V=VOld;
                Theta=ThetaOld;
                Friction=FrictionOld;
                Disp=DispOld;
                CC=0;
            else
                V=VOld;
                Disp=DispOld+(VOld+V)/2*Dt; % Update disp
                % Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; % Ruina Law
                Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; % Dieterich Law
                Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
                Instability=1;
                IsInstability=1;
                CC=0;
            end
        end
    else % Update last
        % =================================================================
        Disp=DispOld+(VOld+V)/2*Dt;
        % Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; % Ruina Law
        Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; % Dieterich Law
        Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
        CC=abs((VTest./V)-1); % Update Critirion
    end

    if rem(Iteration,500)==0
        CCrit0=CCrit0*2;
    end
end
Return_V=V;
Return_Friction=Friction;
Return_Disp=Disp;
Return_Theta=Theta;
Return_ESigma=ESigma;
Return_Dt=Dt;
Return_InstabilityThistime=IsInstability;
% if Disp-Disp1>1e-6
%     fprintf('# A rough calculation\n');
% end
%
% for kk=1:length(V)
%     if V(kk)<1e-20;
%         V(kk)=1e-20;
%         Disp(kk)=DispOld(kk);
%         Theta(kk)=ThetaOld(kk);
%         Friction(kk)=FrictionOld(kk);
%         Dt=Dt/5;
%     end
% end
