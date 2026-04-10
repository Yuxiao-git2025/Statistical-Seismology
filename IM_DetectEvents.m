% Here, all the segments on the same fault are summed up for seismic moment
% and it is judged whether the total rate is greater than the threshold
% rate. If so, it can be determined whether there is an earthquake on this
% fault, rather than checking for the occurrence of an earthquake on each
% segment
function [Event,Fault,CumEvent]=IM_DetectEvents(ITime,IDisp)
% Outputs
%   Event:
%       .Count, .Time, .Location, .Moment, .Mw, .FaultID, .BeginIdx, .EndIdx
%   Fault:
%       .MomentAccum [Fseg x Step], .MomentRate [Fseg x Step]
%   CumEvent: [Step x 1] event-count
global Fseg Fcen flen fID fcount G MomentCrit MagCrit Width
clear Event Fault CumEvent
Step=length(ITime);

% >> Fault moment accumulation and rate
MomentCum=zeros(Fseg,Step);
for iseg=1:fcount
    ID=fID(iseg);
    % Moment release: M0=G*(L*CritLen)*D
    MomentCum(ID,:)=MomentCum(ID,:)+(IDisp(:,iseg)'*(flen(iseg)*Width*G));
end

MomentRate=zeros(Fseg,Step);
dt=diff(ITime(:)).'; 
for ID=1:Fseg
    MomentRate(ID,2:Step)=diff(MomentCum(ID,:))./dt;
end

% >> Event detection for each fault (thresholding on moment rate)
EventTime=[];
EventLocation=[];
EventMoment=[];
EventMw=[];
EventFaultID=[];
EventBeginIdx=[];
EventEndIdx=[];

Count=0;
for ID=1:Fseg
    inEvent=false;
    beginIdx=NaN;

    for j=1:Step
        if MomentRate(ID,j)>MomentCrit
            if ~inEvent
                beginIdx=j;
                inEvent=true;
            end
        else
            if inEvent
                endIdx=j;
                inEvent=false;
                M0=MomentCum(ID,endIdx)-MomentCum(ID,beginIdx);
                if M0>0
                    Mw=(log10(M0)-9.05)/1.5;
                    % ---- Accept only if Mw>=MagCrit
                    if Mw>=MagCrit
                        Count=Count+1;
                        EventBeginIdx(Count,1)=beginIdx;
                        EventEndIdx(Count,1)=endIdx;
                        EventTime(Count,1)=ITime(beginIdx);
                        EventLocation(Count,:)=Fcen(ID,:);
                        EventMoment(Count,1)=M0;
                        EventMw(Count,1)=Mw;
                        EventFaultID(Count,1)=ID;
                    end
                end
            end
        end
    end
end

% >> Event count vs time
CumEvent=zeros(Step,1);
if Count>0
    % Faster than double loop: sort times and do running count
    tEvents=sort(EventTime(:));
    k=1;
    for i=1:Step
        while k<=length(tEvents) && tEvents(k)<=ITime(i)
            k=k+1;
        end
        CumEvent(i)=k-1;
    end
end

% output
Event = struct();
Event.Count = Count;
Event.Time = EventTime;
Event.Location = EventLocation;
Event.Moment = EventMoment;
Event.Mw = EventMw;
Event.FaultID = EventFaultID;
Event.BeginIdx = EventBeginIdx;
Event.EndIdx = EventEndIdx;

Fault = struct();
Fault.MomentAccum = MomentCum;
Fault.MomentRate = MomentRate;
fprintf('# Total %d events (>=%.2f)\n',Event.Count,MagCrit);

end
