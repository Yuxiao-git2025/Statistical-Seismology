function [fcenter,flen,fangle,frad,fRL,finiNormal,finiShear,VlStresscopy,fcount,fseg,flensum,fID]=...
    IM_DiscrestFault()
% Notice the dimension of VlStress will extend here, and outside is input
global Fseg Fcen Flen Fangle FRL FiniNormal FiniShear VlStress CritLen
% =========================================================================
Frad=deg2rad(Fangle);
% Find the end points of all faults
Fx1=Fcen(:,1)+Flen(:)/2.*sin(Frad(:)); % right points
Fy1=Fcen(:,2)+Flen(:)/2.*cos(Frad(:)); % up points
% FaultX2_Bulk=Fcen(:,1)-Flen(:)/2.*sin(Frad(:));
% FaultY2_Bulk=Fcen(:,2)-Flen(:)/2.*cos(Frad(:));

fseg=zeros(Fseg,1);
for i=1:Fseg
    if Flen(i) > CritLen
        fseg(i)=ceil(Flen(i)/CritLen);
    else
        fseg(i)=1;
    end
end

id=0;
for i=1:Fseg
    SegIndex=0;
    for j=1:fseg(i)
        id=id+1;
        SegIndex=SegIndex+1;
        if fseg(i)==1
            fcenter(id,1)=Fcen(i,1);
            fcenter(id,2)=Fcen(i,2);
            flen(id,1)=Flen(i);
            flensum(id,1)=Flen(i);
            fangle(id,1)=Fangle(i);
            frad(id,1)=Frad(i);
            fRL(id,1)=FRL(i);
            finiNormal(id,1)=FiniNormal(i);
            finiShear(id,1)=FiniShear(i);
            VlStresscopy(id,1)=VlStress(i,1);
            fID(id,1)=i;
        else
            flen(id,1)=Flen(i)/fseg(i);
            fangle(id,1)=Fangle(i);
            frad(id,1)=Frad(i);

            fx1(id,1)=Fx1(i)-(SegIndex-1)*flen(id).*sin(frad(id));
            fy1(id,1)=Fy1(i)-(SegIndex-1)*flen(id).*cos(frad(id));
            fx2(id,1)=Fx1(i)-SegIndex*flen(id).*sin(frad(id));
            fy2(id,1)=Fy1(i)-SegIndex*flen(id).*cos(frad(id));

            fcenter(id,1)=(fx1(id,1)+fx2(id,1))/2;
            fcenter(id,2)=(fy1(id,1)+fy2(id,1))/2;
            flensum(id,1)=Flen(i);

            fRL(id,1)=FRL(i);
            finiNormal(id,1)=FiniNormal(i);
            finiShear(id,1)=FiniShear(i);
            VlStresscopy(id,1)=VlStress(i);
            fID(id,1)=i;

        end

    end
end
fcount=id;

end

