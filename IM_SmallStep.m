function [Return_V,Return_Friction,Return_Disp,Return_Theta,Return_ESigma,Return_Dt,Return_InstabilityThistime]...
    = IM_SmallStep(CCrit,DispOld,FrictionOld,ThetaOld,VOld,...
    Friction0,a,b,Dc,V0,kShearii,Dt,Omega,G,Instability,Disp, Friction, V, ...
    Theta, DV, Esigma,TotalDisp)
% Here refer to (Im 2017) using the constant friction solution
% Because it needs smaller time-step like dt=1e-5
Iteration=0;
IsInstability=0; % 0/1/3
CC=0.1;
CCrit0=CCrit;

% =========================================================================
while CC>CCrit0
    Iteration=Iteration+1;
    VOldIter=V;
    VTest=V; 

    %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
    Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
    Friction=Friction0+a*log(V/V0)+b*log(Theta*V0./Dc)+G/6000/Esigma*V; % Initial friction
    F=TotalDisp-Friction*Esigma/kShearii;
    Disp=(DispOld-F)*cos(Omega*Dt)+(VOld/Omega)*sin(Omega*Dt)+F; % Equation (9)
    V=(Disp-DispOld)/Dt*2-VOld;

    Fv1=VOldIter-V; % We are testing this Newton Rhapson Function. Lets send this to zero

    % Deviated value for NR
    V=VOldIter+DV;
    VOldIter=V;
    %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
    Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
    Friction=Friction0+a*log(V/V0)+b*log(Theta*V0./Dc)+G/6000/Esigma*V; % Initial friction
    F=TotalDisp-Friction*Esigma/kShearii;
    Disp=(DispOld-F)*cos(Omega*Dt)+(VOld/Omega)*sin(Omega*Dt)+F; % Equation (9)
    V=(Disp-DispOld)/Dt*2-VOld;

    Fv2=VOldIter-V; % We are testing this Newton Rhapson Function. Lets send this to zero
    DF=(Fv2-Fv1)/DV;
    V=VTest-Fv1/DF; % Update velocity

    % =====================================================================
    if V<0
        if Instability==1
            V=1e-20;
            Disp=DispOld+(VOld+V)/2*Dt; % Update disp
            %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
            Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
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
                %                 V=1e-10;
                %                 V=VOld;
                if V>1e-20
                    V=VOld/2;
                else
                    V=1e-20;
                    Instability=1;
                    IsInstability=1;
                end
                Disp=DispOld+(VOld+V)/2*Dt; % Update disp
                %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
                Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
                Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
                CC=0;
            end
        end
    elseif ~isreal(V)
        if Instability==1
            V=1e-20;
            Disp=DispOld+(VOld+V)/2*Dt; % Update disp
            %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
            Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
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
                V=VOld;
                Disp=DispOld+(VOld+V)/2*Dt; % Update disp
                %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
                Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
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
            %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
            Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
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
                %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
                Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
                Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
                Instability=1;
                IsInstability=1;
                CC=0;
            end
        end
    else
        % =================================================================
        Disp=DispOld+(VOld+V)/2*Dt; % Update disp
        %         Theta=ThetaOld-V*ThetaOld/Dc*log(V*ThetaOld/Dc)*Dt; %Ruina Law
        Theta=ThetaOld+(1-V*ThetaOld/Dc)*Dt; %Dieterich Law
        Friction=Friction0+a.*log(V./V0)+b.*log(Theta*V0./Dc);
        CC=abs(VTest./V-1);
    end

    if rem(Iteration,100)==0
        CCrit0=CCrit0*10;
    end
    if CCrit0>1e-7
        IsInstability=2;
    end
end
Return_V=V;
Return_Friction=Friction;
Return_Disp=Disp;
Return_Theta=Theta;
Return_ESigma=Esigma;
Return_Dt=Dt;
Return_InstabilityThistime= IsInstability;

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
