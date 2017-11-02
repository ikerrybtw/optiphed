function [success,dmHE,dmIPOX,SUFF]=calculateDistanceOfHandENucleiToIPOX(dmHE, dmIPOX,xI,yI, nCells, ipoxIntensity)
%%heI:=random subset of H&E nuclei (for performance reasons)
% rng(2)
global MAXSECONDBESTPEAK
MAXSECONDBESTPEAK=50;%85;%65;

MAXDIST=0.0075; %0.05
SUFF='_ClosestIPOX';
if ~exist('nCells','var')
    heI=1:size(dmHE,1);
else
    %%Find IPOX nucleus of high staining absorption
    q=quantile(double(dmIPOX(:,'Intensity_MaxIntensityEdge_Costum')),ipoxIntensity);
    [~,ia]=min(abs(dmIPOX(:,'Intensity_MaxIntensityEdge_Costum')-q));
    x=pdist2(double(dmIPOX(ia,{xI,yI})),double(dmHE(:,{xI,yI})));
    [~,ia]=sort(x);
    %%Pick nCells closest to IPOX nucleus of high staining absorption
    center=double(dmHE(ia(1),{xI,yI}));
    x=pdist2(center,double(dmHE(:,{xI,yI})));
    [~,ia]=sort(x);
    heI=ia(1:nCells);
    %%Extent to form square:
    ia=find(dmHE(:,xI)>=min(dmHE(heI,xI))-MAXDIST & dmHE(:,xI)<=max(dmHE(heI,xI))+MAXDIST);
    ib=find(dmHE(:,yI)>=min(dmHE(heI,yI))-MAXDIST & dmHE(:,yI)<=max(dmHE(heI,yI))+MAXDIST);
    ij=intersect(ia,ib);
    heI=union(heI,ij);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%FINE TUNING ALIGNMENT H&E vs. IPOX%%%%%%%%

%Fine tuning
try
    
    [dmHE,dmIPOX,success]=callFineTuneAlignment(dmHE,heI,dmIPOX,1:0.5:2,MAXDIST,xI,yI, ipoxIntensity);
    if ~success
        return;
    end
catch err
    disp(err)
    success=0;
    return;
end

% %%Limit H&E and IPOX XY-coordinates
% ia=find(dmIPOX(:,xI)>=0.5-MAXDIST & dmIPOX(:,xI)<=0.5+MAXDIST);
% ib=find(dmIPOX(:,yI)>=0.5-MAXDIST & dmIPOX(:,yI)<=0.5+MAXDIST);
% dmIPOX=dmIPOX(intersect(ia,ib),:);
% ia=find(dmHE(:,xI)>=0.5-MAXDIST & dmHE(:,xI)<=0.5+MAXDIST);
% ib=find(dmHE(:,yI)>=0.5-MAXDIST & dmHE(:,yI)<=0.5+MAXDIST);
% dmHE=dmHE(intersect(ia,ib),:);

% OLDCOLS=unique({'AreaShape_Area','AreaShape_Compactness', 'AreaShape_Eccentricity', 'AreaShape_Extent', 'AreaShape_FormFactor', 'AreaShape_MajorAxisLength', 'AreaShape_MaxFeretDiameter', 'AreaShape_MaximumRadius', 'AreaShape_MeanRadius', 'AreaShape_MinFeretDiameter', 'AreaShape_MinorAxisLength', 'AreaShape_Orientation', 'AreaShape_Perimeter', 'AreaShape_Zernike_0_0','AreaShape_Zernike_1_1', 'AreaShape_Zernike_3_3',xI,yI});
NEWCOLS={};%strcat(OLDCOLS,SUFF);
NEWCOLS=setdiff([NEWCOLS,'ClosestIPOXNuclei','ClosestHandENuclei','N_Neighbors_1E-3'],colnames(dmHE));
for col = NEWCOLS
    dmHE=addColumn(dmHE,char(col),nan);
    dmHE(:,col)=nan;
end

% dmHE_=dmHE(heI,:);
for i=1:size(dmHE,1)
    if mod(i,500)==0
        disp(['Processed ',num2str(i),' out of ',num2str(size(dmHE,1)),' nuclei...']);
    end
    x=pdist2(double(dmHE(i,{xI,yI})),double(dmIPOX(:,{xI,yI})));
    [a,ia]=nanmin(x);
    if isempty(ia)
        continue;
    end
    dmHE(i,'ClosestIPOXNuclei')=a;
    %     dmHE(i,strcat(OLDCOLS,SUFF))=dmIPOX_(ia,OLDCOLS);
    dmHE(i,'Intensity_MaxIntensityEdge_Costum')=dmIPOX(ia,'Intensity_MaxIntensityEdge_Costum');
    dmHE(i,'Intensity_UpperQuartileIntensity_Costum')=dmIPOX(ia,'Intensity_UpperQuartileIntensity_Costum');
    y=pdist2(double(dmHE(i,{xI,yI})),double(dmHE(setdiff(1:size(dmHE,1),i),{xI,yI})));
    dmHE(i,'ClosestHandENuclei')=nanmin(y);
    dmHE(i,'N_Neighbors_1E-3')=length(find(y<MAXDIST*2));
    %     disp(dmHE(i,:))
end
% rng('default')
end


function [dmHE,dmIPOX,success]=callFineTuneAlignment(dmHE,heI,dmIPOX,zooms,MAXDIST,xI,yI, ipoxIntensity)
global MAXSECONDBESTPEAK
import bioma.data.*
success=1;
xyAdjust=[];percentDecline=[];

xoyo=mean(dmIPOX(:,{xI,yI})); 
dmIPOXRot=dmIPOX;
for a=[strcat({'clockwise '},cellstr(num2str([0:0.5:3]'))); strcat({'anticlockwise '},cellstr(num2str([0.5:0.5:3]')))]' %%Try slightly different angles
    tmp=strsplit(char(a),' ');
    theta=str2num(tmp{2});   rotDirection=tmp{1};
    [xr,yr] = rotateData(double(dmIPOX(:,xI)),double(dmIPOX(:,yI)),...
        xoyo(1),xoyo(2),degtorad(theta),char(rotDirection));
    dmIPOXRot(:,{xI,yI})=[xr;yr]';
    
    
    for m=zooms %%Zoom into variable scales
        %Limit H&E XY-coordinates
        ia=find(dmHE(:,xI)>=min(dmHE(heI,xI))-MAXDIST*m & dmHE(:,xI)<=max(dmHE(heI,xI))+MAXDIST*m);
        ib=find(dmHE(:,yI)>=min(dmHE(heI,yI))-MAXDIST*m & dmHE(:,yI)<=max(dmHE(heI,yI))+MAXDIST*m);
        ij=intersect(ia,ib);
        %Limit IPOX XY-coordinates
        ia=find(dmIPOXRot(:,xI)>=min(dmHE(ij,xI))-MAXDIST*2*m & dmIPOXRot(:,xI)<=max(dmHE(ij,xI))+MAXDIST*2*m);
        ib=find(dmIPOXRot(:,yI)>=min(dmHE(ij,yI))-MAXDIST*2*m & dmIPOXRot(:,yI)<=max(dmHE(ij,yI))+MAXDIST*2*m);
        ipoxI=intersect(ia,ib);
        
        %Check if these have been realigned already:
        testCheck=dmIPOX(ipoxI,{[xI,'_NoLocalAdjustment'],[yI,'_NoLocalAdjustment']})==dmIPOX(ipoxI,{xI,yI});
        if ~all(all(testCheck))
            disp('A subset of these coordinates did already undergo local realignment. Skipping');
            success=0;
            return;
        end
        
        [xyAdjust_,percentDecline_, dmHE1, dmIPOX1]=finetuneAlignment(dmHE(ij,:),dmIPOXRot(ipoxI,:),xI,yI, ipoxIntensity);
%         if percentDecline_<MAXSECONDBESTPEAK
%             dmHE=dmHE1; dmIPOX=dmIPOX1;
%             return;
%         end
        xyAdjust=[xyAdjust; xyAdjust_];
        percentDecline=[percentDecline,percentDecline_];
    end
    
end
%     %%Adjust IPOX coordinates
%     disp(['Final X-offset: ',num2str(mean(xyAdjust(:,1)))])
%     disp(['Final Y-offset: ',num2str(mean(xyAdjust(:,2)))])
%     disp(['Final angle: ',char(finalA)])
%     dmIPOX_(:,xI)=dmIPOX_(:,xI)+mean(xyAdjust(:,1));
%     dmIPOX_(:,yI)=dmIPOX_(:,yI)+mean(xyAdjust(:,2));


if all(percentDecline>=MAXSECONDBESTPEAK)
    disp(['Matching IPOX to H&E nuclei was not succesfull (%decline=',num2str(percentDecline),'). Skipping assignment.'])
    success=0;
    return;
elseif length(zoom)==1
    dmHE=dmHE1; dmIPOX=dmIPOX1;
else
    [dmHE,dmIPOX,success]=callFineTuneAlignment(dmHE,heI,dmIPOX,zooms(percentDecline==min(percentDecline)),MAXDIST,xI,yI, ipoxIntensity);
end
end





