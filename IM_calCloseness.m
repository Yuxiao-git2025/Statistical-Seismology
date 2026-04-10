% This function removes the faults with too strong interaction 
function [fcen1,flen1,fangle1,fRL1,finiNormal1,finiShear1,fcount1,flensum1,fID1,VlStress1,fseg1]...
    =IM_calCloseness(kShear,kNormal)
global fcen flen fangle fRL finiNormal finiShear fcount flensum fID VlStress
global Doing
TooClose=zeros(fcount,1); % 0 for not, 1 for yes
UnstableCount=0;
EShear=zeros(fcount,fcount);
ENormal=zeros(fcount,fcount);
for i=1:fcount
    for j=1:fcount
        if i==j
            EShear(i,j)=0; % do not consider interaction itself
        else
            EShear(i,j)=kShear(i,j)/kShear(i,i);
            ENormal(i,j)=kNormal(i,j)/kShear(i,i);
        end
    end
end


for i=1:fcount
    for j=1:fcount
        if abs(EShear(i,j))>0.8
            if TooClose(i)==0 && TooClose(j)==0
                if flensum(i)<flensum(j)
                    UnstableCount=UnstableCount+1;
                    TooClose(i)=1;   
                else
                    UnstableCount=UnstableCount+1;
                    TooClose(j)=1;
                end
            end
        end
        if abs(ENormal(i,j))>0.8
            if TooClose(i)==0 && TooClose(j)==0
                if flensum(i)<flensum(j)
                    UnstableCount=UnstableCount+1;
                    TooClose(i)=1;
                else
                    UnstableCount=UnstableCount+1;
                    TooClose(j)=1;
                end
            end
        end
    end
end

fcount0=0;
for i=1:length(TooClose)
    if TooClose(i)==0
        fcount0=fcount0+1;
        fcen0(fcount0,:)=fcen(i,:);
        flen0(fcount0,1)=flen(i);
        fangle0(fcount0,1)=fangle(i);
        fRL0(fcount0,1)=fRL(i);
        finiShear0(fcount0,1)=finiShear(i);
        finiNormal0(fcount0,1)=finiNormal(i);
        flensum0(fcount0,1)=flensum(i);
        fID0(fcount0,1)=fID(i);
        VlStress0(fcount0,1)=VlStress(i,1);
    end
end
fcount1=fcount0;
fcen1=fcen0;
flen1=flen0;
fangle1=fangle0;
fRL1=fRL0;
finiShear1=finiShear0;
finiNormal1=finiNormal0;
flensum1=flensum0;
fID1=fID0;
VlStress1=VlStress0;
[~,~,fidCompact]=unique(fID1);
fseg1=accumarray(fidCompact,1);
% Fseg1=length(fseg1); % do not update ! (cause' initial fault number)

% Update main/background counts after removing too-close elements
global FnumMain fnumMain FnumBkgd fnumBkgd
if Doing==2
    global Fseg MainPart
    nMain=round(Fseg*MainPart);
    % (1) element-level counts
    fnumMain=sum(fID1<=nMain);
    fnumBkgd=sum(fID1>nMain);
    % (2) fault-level counts
    ID1=unique(fID1);
    FnumMain=sum(ID1<=nMain);    % remaining main-band faults
    FnumBkgd=sum(ID1>nMain);    % remaining background faults
else 
    % Here only consider the main-faults
    FnumMain=length(fseg1);
    FnumBkgd=0;
    fnumMain=fcount1; % or sum(fseg1)
    fnumBkgd=0;
end
end

