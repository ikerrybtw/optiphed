function [xyAdjust,kurto,dmHE_,dmIPOX_]=finetuneAlignment(dmHE,dmIPOX,xI,yI,ipoxIntensity)
global MAXSECONDBESTPEAK
%%Convert to images
MSIZE=40;
if size(dmHE,1)>1000
    MSIZE=ceil(2000/(log(size(dmHE,1)+500)^power(log(size(dmHE,1)+500),1/3)))-1;
end
%MSIZE=ceil(60/(log10(size(dmHE,1))^2)); %12/log10(size(dmHE,1))); 
dmIPOX=removeNoiseFromImage(dmIPOX,xI,yI);
% feature='AreaShape_MajorAxisLength';
f1=figure('name',['IPOX',strrep(num2str(ipoxIntensity),'  ','-')],'units','normalized','Position',[0,0,2.25,1.5]); %
plot(double(dmIPOX(:,xI)),double(dmIPOX(:,yI)),'.','MarkerSize',MSIZE,'MarkerEdgeColor','k','MarkerFaceColor','k');
set(gca,'XLim',[nanmin(dmIPOX(:,xI)),nanmax(dmIPOX(:,xI))],'YLim',[nanmin(dmIPOX(:,yI)),nanmax(dmIPOX(:,yI))]);
axis off; hold on;
% addColors(dmIPOX,feature,xI,yI,MSIZE);
ipoxF = getframe(gca);
 
xL=get(gca,'XLim'); yL=get(gca,'YLim');
ipoxImgCenter=[mean(xL),mean(yL)];
f2=figure('name',['HandE_',strrep(num2str(ipoxIntensity),'  ','-')],'units','normalized','Position',[0,0,2.25,1.5]); 
plot(double(dmHE(:,xI)), double(dmHE(:,yI)),'.','MarkerSize',MSIZE,'MarkerEdgeColor','k','MarkerFaceColor','k'); set(gca,'XLim',xL,'Ylim', yL); axis off; hold on;
% addColors(dmHE,feature,xI,yI,MSIZE);
heF = getframe(gca);

%Calculate offset
slides=struct('ipox',ipoxF.cdata(:,:,2)+ipoxF.cdata(:,:,3),'hande',heF.cdata(:,:,2)+heF.cdata(:,:,3));
for f=fieldnames(slides)'
    template=getfield(slides,char(f));
    xD=find(var(double(template),[],2)==0);  template(xD,:)=[];
    yD=find(var(double(template))==0);       template(:,yD)=[];
    slides=setfield(slides,char(f),template);
end

crr = xcorr2(slides.ipox,slides.hande);
[ssr,snd] = max(crr(:)); crrsub=crr(round(snd*0.85):round(snd*1.15));
[ij,ji] = ind2sub(size(crr),snd); %%the lower-right corner of the IPOX section that corresponds to the lower right corner of the H&E section
f5=figure('name',['CrossCorrMax_',strrep(num2str(ipoxIntensity),'  ','-')]);
kurto=100*max(crrsub([1:1000,(length(crrsub)-1000):length(crrsub)]))/max(crrsub);
plot(crr(:)); title(['MinAsPercentOfMax: ',num2str(kurto,2)]); hold on
plot(snd,ssr,'or'); plot([round(snd*0.85),round(snd*1.15)],[ssr,ssr],'b*'); hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Save normalized matching coordinates%%%%
dmIPOX_=dmIPOX; dmHE_=dmHE;
dmIPOX_(:,{xI,yI})=normalize01(double(dmIPOX(:,{xI,yI})));
dmHE_(:,{xI,yI})=normalize01(double(dmHE(:,{xI,yI})));

% %%MAP back
xyRange=[max(dmIPOX(:,xI))-min(dmIPOX(:,xI)),max(dmIPOX(:,yI))-min(dmIPOX(:,yI)) ];
% yRange=[ij-size(slides.hande,1)+1, ij]./size(slides.ipox,1);%% The % of the range
% yRange=max(dmIPOX(:,yI))-yRange*xyRange(2);
% xRange=[ji-size(slides.hande,2)+1, ji]./size(slides.ipox,2);
% xRange=max(dmIPOX(:,xI))-xRange*xyRange(1);
% f4=figure('name',['mapBack_',strrep(num2str(ipoxIntensity),'  ','-')]);
% ia=find(dmIPOX(:,xI)>=xRange(2) & dmIPOX(:,xI)<=xRange(1));
% ib=find(dmIPOX(:,yI)>=yRange(2) & dmIPOX(:,yI)<=yRange(1));
% iEq=intersect(ia,ib);
% plot(double(dmIPOX(iEq,xI)),double(dmIPOX(iEq,yI)),'.','MarkerSize',MSIZE,'MarkerEdgeColor','k','MarkerFaceColor','k'); axis off;


%%TO map back: 1) convert impox image into 2D coordinates, while marking the
%%H&E-mapping section as 3rd column; 2) normalize coordinates between 0 and
%%1: H&E, IPOX and converted IPOX. 3) Finally select only those IPOX nuclei
%%within the specified matching rectangle
MAPMARKERID=9999;
tmp=double(slides.ipox); tmp(ij:-1:ij-size(slides.hande,1)+1,ji:-1:ji-size(slides.hande,2)+1)=MAPMARKERID; %1
ipox_xyCoord=img2Cart(tmp); %1 
ipox_xyCoord(:,1)=normalize01(ipox_xyCoord(:,1)); ipox_xyCoord(:,2)=normalize01(ipox_xyCoord(:,2)); %2
ipox_xyCoord=ipox_xyCoord(ipox_xyCoord(:,3)==MAPMARKERID,:); %3
iEq=find(dmIPOX_(:,xI)>=min(ipox_xyCoord(:,1)) & dmIPOX_(:,xI)<=max(ipox_xyCoord(:,1)) & ...
    dmIPOX_(:,yI)>=min(ipox_xyCoord(:,2)) & dmIPOX_(:,yI)<=max(ipox_xyCoord(:,2)));
f4=figure('name',['mapBack_',strrep(num2str(ipoxIntensity),'  ','-')],'units','normalized','Position',[0,0,0.75,0.75]);
plot(double(dmIPOX(iEq,xI)),double(dmIPOX(iEq,yI)),'.','MarkerSize',MSIZE,'MarkerEdgeColor','k','MarkerFaceColor','k'); axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Calculateing the actual shift%%%%%%
% % %% Alternative 1: Find H&E nucleus closest to the lower-left-most IPOX nucleus (<-- because IPOX+ nuclei are much less uniformly distributed than H&E)
% lrIpox=min(dmIPOX(iEq,{xI,yI}));
% [~,ia]=min(pdist2(lrIpox,double(dmHE(:,{xI,yI}))));
% %%Distance between that H&E nucleus and lower-left-most H&E nucleus represents the shift distance
% xyAdjust=double(dmHE(ia,{xI,yI})-min(dmHE(:,{xI,yI})));
% % Alternative 2:
qs=[0:0.0005:1]; 
% xyAdjust=quantile(dmHE(:,{xI,yI}),qs)-quantile(dmIPOX(iEq,{xI,yI}),qs);
% xyAdjust=median(xyAdjust);
% %% Alternative 3: centering ? expecting center of entire IPOX to also be center of IPOX subset
% xLsub=quantile(dmIPOX(iEq,xI),[0,1]);
% yLsub=quantile(dmIPOX(iEq,yI),[0,1]);
% subIpoxImgCenter=[mean(xLsub),mean(yLsub)];
% xyAdjust=subIpoxImgCenter-ipoxImgCenter;
% %% Alternative 4: - but change in direction seems to be one-sided - how come it works?
% xyAdjust=[ji-size(slides.hande,2)+1,ij-size(slides.hande,1)+1];
% %%Alternative 5 (similar to Alternative 2: 
% subIpoxImgCenter=mean(dmIPOX(iEq,{xI,yI}));
% heImgCenter=mean(dmHE(:,{xI,yI}));
% xyAdjust=heImgCenter-subIpoxImgCenter;
% %%Alternative 6
tmp=double(dmIPOX_(iEq,{xI,yI}))-normalize01(double(dmIPOX_(iEq,{xI,yI})));
xyAdjust=-1*mean(tmp);
xyAdjust=xyAdjust.*xyRange;

slides.ipox(ij:-1:ij-size(slides.hande,1)+1,ji:-1:ji-size(slides.hande,2)+1) = rot90(slides.hande,2);
f3=figure('name',['HandE_On_IPOXbackground_IPOXQ_',strrep(num2str(ipoxIntensity),'  ','-')]);
imshow(slides.ipox);

% if kurto <MAXSECONDBESTPEAK
    savefig(f3,[get(f3,'name'),'.fig'],'compact');
    savefig(f5,[get(f5,'name'),'.fig'],'compact');
% end

%%%%%% Save normalized matching coordinates%
dmIPOX_=dmIPOX_(iEq,:);
dmIPOX_(:,{xI,yI})=normalize01(double(dmIPOX_(:,{xI,yI})));
 
disp(['X-offset: ',num2str(xyAdjust(1))])
disp(['Y-offset: ',num2str(xyAdjust(2))]) 
close(f1); close(f2); close(f3); close(f4); close(f5); 
disp('Done Fine tuning')   
end

function addColors(dmIPOX,feature,xI,yI,MSIZE)
cmap=hot(1+round(max(dmIPOX(:,feature))-min(dmIPOX(:,feature))));
allC=1+double(round(dmIPOX(:,feature)-min(dmIPOX(:,feature))));
for c=unique(allC)'
    i=find(c==allC);
    plot(double(dmIPOX(i,xI)),double(dmIPOX(i,yI)),'.','MarkerSize',MSIZE,'MarkerEdgeColor',cmap(c,:),'MarkerFaceColor',cmap(c,:));
end
end



