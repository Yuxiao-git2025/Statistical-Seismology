function [kShear,kNormal]=IM_calStiffness(fcount,fcen,flen,fangle,G,nu,fRL)
% 1. kNormal: The increment of the "normal stress" on segment i 
%       caused by the unit slip of segment j;
% 2. kNormal: The increment of the "shear stress" on segment i 
%       caused by the unit slip of segment j;
% 3. kShearii (Self-stiffness): The influence of self-slip of segment i on 
%       its own shear stress, representing the diagonal elements;
kShear=zeros(fcount,fcount);
kNormal=zeros(fcount,fcount);
% Total elements is nxn; For each element i, the effect of all j elements
% on it must be calculated
for i=1:fcount
    for j=1:fcount
        Thecenter=fcen(j,:);
        Thelen=flen(j)/2;
        Theangle=fangle(j); % Degree from North

        RecCenter=fcen(i,:);
        RecAngle=fangle(i);

        TheDisp=-1;

        nFaultAngle=Theangle+90;
        lFaultAngle=Theangle;
        nFault=cos(deg2rad(nFaultAngle));
        lFault=cos(deg2rad(lFaultAngle));

        nReceiverAngle=RecAngle+90;
        lReceiverAngle=RecAngle;
        nReceiver=cos(deg2rad(nReceiverAngle));
        lReceiver=cos(deg2rad(lReceiverAngle));


        Xi=nFault*(RecCenter(1)-Thecenter(1))-lFault*(RecCenter(2)-Thecenter(2));   % twocurve X
        Zeta=lFault*(RecCenter(1)-Thecenter(1))+nFault*(RecCenter(2)-Thecenter(2)); % onecurve Z

        F3=((nFault^2-lFault^2)*Zeta-2*nFault*lFault*(Xi+Thelen))/((Xi+Thelen)^2+Zeta^2) ...
            - ((nFault^2-lFault^2)*Zeta-2*nFault*lFault*(Xi-Thelen))/((Xi-Thelen)^2+Zeta^2);
        F4=-(2*nFault*lFault*Zeta+(nFault^2-lFault^2)*(Xi+Thelen))/((Xi+Thelen)^2+Zeta^2) ...
            +(2*nFault*lFault*Zeta+(nFault^2-lFault^2)*(Xi-Thelen))/((Xi-Thelen)^2+Zeta^2);
        F5=(nFault*(nFault^2-3*lFault^2)*((Xi+Thelen)^2-Zeta^2)+ ...
            2*lFault*(3*nFault^2-lFault^2)*(Xi+Thelen)*Zeta)/...
            ((Xi+Thelen)^2+Zeta^2)^2 ...
            -(nFault*(nFault^2-3*lFault^2)*((Xi-Thelen)^2-Zeta^2)+ ...
            2*lFault*(3*nFault^2-lFault^2)*(Xi-Thelen)*Zeta)/...
            ((Xi-Thelen)^2+Zeta^2)^2;
        F6=(2*nFault*(nFault^2-3*lFault^2)*(Xi+Thelen)*Zeta-lFault ...
            *(3*nFault^2-lFault^2)*((Xi+Thelen)^2-Zeta^2))/((Xi+Thelen)^2+Zeta^2)^2 ...
            -(2*nFault*(nFault^2-3*lFault^2)*(Xi-Thelen)*Zeta-lFault ...
            *(3*nFault^2-lFault^2)*((Xi-Thelen)^2-Zeta^2))/((Xi-Thelen)^2+Zeta^2)^2;

        StressXX=TheDisp*G/2/pi/(1-nu)*(2*nFault^2*F3 -2*nFault*lFault*F4 +Zeta*(nFault*F5 -lFault*F6));
        StressZZ=-TheDisp*G/2/pi/(1-nu)*(2*lFault^2*F3 +2*nFault*lFault*F4 + Zeta*(nFault*F5 -lFault*F6));
        StressXZ=(TheDisp*G/2/pi/(1-nu)*(F4+Zeta*(lFault*F5+nFault*F6)));

        if abs(StressXX)>1e10; StressXX=0; end
        if abs(StressZZ)>1e10; StressZZ=0; end
        if abs(StressXZ)>1e10; StressXZ=0; end
        % The shear-stiffness dependent not only receiver i but impactings
        kShear(i,j)=fRL(i)*fRL(j)*(nReceiver*lReceiver*(StressXX-StressZZ)+(nReceiver^2-lReceiver^2)*StressXZ);
        % Tension is positive and compression is negative
        kNormal(i,j)=fRL(j)*(lReceiver^2*StressXX+nReceiver^2*StressZZ+2*nReceiver*lReceiver*StressXZ);

        %         ResultSXX(i,j)=StressXX;
        %         ResultSZZ(i,j)=StressZZ;
        %         ResultSXZ(i,j)=StressXZ;
        %         ResultReceiverShear(i,j)=ReceiverRRLL*(nReceiver*lReceiver*
        %               (StressXX-StressZZ)+(nReceiver^2-lReceiver^2)*StressXZ);
        %         ResultReceiverNormal(i,j)=-(lReceiver^2*StressXX+nReceiver^2
        %               *StressZZ+2*nReceiver*lReceiver*StressXZ);
    end
end

end

