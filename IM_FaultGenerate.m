function [fcen, flen, fangle, frad, fRL, finiNormal, finiShear, fcount, ...
    fseg, flensum, fID, Fcen, Flen, Fangle, Frad, FRL, FiniNormal, FiniShear] ...
    =IM_FaultGenerate()
global Fseg Fangle1 Fangle2 FlenMin FlenMax DxMin DxMax DyMin DyMax 
global Sigma1 Sigma3 
global MainPart CritLen FnumMain fnumMain FnumBkgd fnumBkgd 
global FnumMainini FnumBkgdini fnumMainini fnumBkgdini
global seed
% =========================================================================
% Input Parameters:
rng(seed);
% Main-fault Range in XAxis 
% (also can be put in Global function, but in here you can specify)
xRange=100;
% Main-fault Range in YAxis
MainYmin=0.2*DyMax; 
MainYmax=0.8*DyMax;


% =========================================================================
% ============= Generate bulk fault lengths (main + BG) =============
RandLen = slicesample(1, Fseg, 'pdf', @exppdf, 'thin', 5, 'burnin', 1e3);
RandLen = RandLen / max(RandLen);
Flen = 10.^( RandLen*(log10(FlenMax)-log10(FlenMin)) + log10(FlenMin) );

RandBG = slicesample(1, Fseg, 'pdf', @exppdf, 'thin', 5, 'burnin', 1e3);
RandBG = RandBG / max(RandBG);
FaultLenBulk_BG = 10.^( RandBG*(log10(FlenMax)-log10(FlenMin)) + log10(FlenMin) );

% ============= Generate bulk fault centers (main + BG) =============
% Main band centers: use repulsive Y then map to [MainYmin, MainYmax]
[CenterRep, Lrep] = IM_Repulsive(Fseg);
FaultCenterY = CenterRep(:,2) * (MainYmax - MainYmin) / Lrep + MainYmin;

Ymid = (MainYmax + MainYmin)/2;
Fcen = zeros(Fseg,2);
Fcen(:,2) = FaultCenterY;
Fcen(:,1) = (FaultCenterY - Ymid) * tan(deg2rad(Fangle1)) + Ymid + ...
                    normrnd(0, xRange, [Fseg,1]);

% Background centers: repulsive in box -> map to [DxMin,DxMax]x[DyMin,DyMax]
[CenterRepBG, LrepBG] = IM_Repulsive(Fseg);
FaultCenBulk_BG = zeros(Fseg,2);
FaultCenBulk_BG(:,1) = CenterRepBG(:,1) * (DxMax - DxMin) / LrepBG + DxMin;
FaultCenBulk_BG(:,2) = CenterRepBG(:,2) * (DyMax - DyMin) / LrepBG + DyMin;

% Replace tail portion with BG
nMain = round(Fseg*MainPart);
FnumMain = nMain;
FnumBkgd = Fseg-nMain;
Fcen(nMain+1:end,:) = FaultCenBulk_BG(nMain+1:end,:);
Flen(nMain+1:end)   = FaultLenBulk_BG(nMain+1:end);

Fcen(:,1) = min(max(Fcen(:,1), DxMin), DxMax);
Fcen(:,2) = min(max(Fcen(:,2), DyMin), DyMax);

% ============= Assign bulk angles (Mix Fangle1/Fangle2) and stresses =============
FRand=2*rand(Fseg,1)-1;
Fangle=zeros(Fseg,1);
Frad=zeros(Fseg,1);
FRL=zeros(Fseg,1);
FiniNormal=zeros(Fseg,1);
FiniShear=zeros(Fseg,1);
% Each fault' angle randomly obtained from Fangle1 and Fangle2
for i = 1:Fseg
    if FRand(i)>0.8
        Fangle(i)=Fangle1;
        FRand(i)=1;
    else
        Fangle(i)=Fangle2;
        FRand(i)=-1;
    end
    Frad(i) = deg2rad(Fangle(i));
    if Fangle(i) > 0
        FRL(i) = -1;
    else
        FRL(i) = 1;
    end
    % Calculate the initial stress
    FiniNormal(i) = (Sigma1+Sigma3)/2 + (Sigma1-Sigma3)/2 * cos(2*(pi/2 - Frad(i)));
    FiniShear(i) = abs((Sigma1-Sigma3)/2 * sin(2*(pi/2 - Frad(i))));
end

% ============= Bulk endpoints (For segmentation) =============
X1 = Fcen(:,1) + Flen/2 .* sin(Frad);
X2 = Fcen(:,1) - Flen/2 .* sin(Frad);
Y1 = Fcen(:,2) + Flen/2 .* cos(Frad);
Y2 = Fcen(:,2) - Flen/2 .* cos(Frad);

% ============= Decide segment =============
fseg = ones(Fseg,1);
for i = 1:Fseg
    if Flen(i) > CritLen
        fseg(i) = ceil(Flen(i)/CritLen);
    end
end

% ============= Discretize into elements =============
fSum = sum(fseg);
fcen   = zeros(fSum,2);
flen   = zeros(fSum,1);
fangle = zeros(fSum,1);
frad   = zeros(fSum,1);
fRL    = zeros(fSum,1);
finiNormal = zeros(fSum,1);
finiShear  = zeros(fSum,1);
flensum = zeros(fSum,1);
fID = zeros(fSum,1);

ID=0;
for i=1:Fseg
    nseg=fseg(i);
    for j=1:nseg
        ID=ID+1;
        fID(ID)=i;
        fangle(ID)=Fangle(i);
        frad(ID)=Frad(i);
        fRL(ID)=FRL(i);
        finiNormal(ID)=FiniNormal(i);
        finiShear(ID)=FiniShear(i);
        flensum(ID)=Flen(i);

        if nseg == 1
            fcen(ID,:) = Fcen(i,:);
            flen(ID)   = Flen(i);
        else
            segLen = Flen(i)/nseg;
            flen(ID) = segLen;

            % Same discretization as your original:
            % start from (X1,Y1) and march towards (X2,Y2)
            x1 = X1(i) - (j-1)*segLen*sin(Frad(i));
            y1 = Y1(i) - (j-1)*segLen*cos(Frad(i));
            x2 = X1(i) - j*segLen*sin(Frad(i));
            y2 = Y1(i) - j*segLen*cos(Frad(i));

            fcen(ID,:) = [(x1+x2)/2, (y1+y2)/2];
        end
    end
end
fcount=ID;
fcen = fcen(1:fcount,:);
flen = flen(1:fcount);
fangle = fangle(1:fcount);
frad = frad(1:fcount);
fRL = fRL(1:fcount);
finiNormal = finiNormal(1:fcount);
finiShear = finiShear(1:fcount);
flensum = flensum(1:fcount);
fID = fID(1:fcount);
% Find the main-band fault-segments numbers
fnumMain = sum(fID<=nMain);
fnumBkgd = sum(fID>nMain);

% Resave the faults number:
FnumMainini=FnumMain;
FnumBkgdini=FnumBkgd;
fnumMainini=fnumMain;
fnumBkgdini=fnumBkgd;
end
