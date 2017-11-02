function slide_dm=overlayAdjacentSlides(slide_dm, base_Image,XCOORDLAB, YCOORDLAB, subI, MAXITER)
%%% returns input matrix with nuclei coordinates adjusted such that cells
%%% match/overlay well between images
%% images are aligned to match first image (i.e. nuclei with smalles ImageNumber entries)

% cellProfilerOutF='1820_T_105452F_Nuclei_Rout.txt';
% base_Image='../1820_T_105452F_1_CDX2_layer4.tif'
import bioma.data.*

if ~exist('MAXITER','var')
    MAXITER=4;
end
if ~exist('XCOORDLAB','var')
    XCOORDLAB='AreaShape_Center_X';
end
if ~exist('YCOORDLAB','var')
    YCOORDLAB='AreaShape_Center_Y';
end
SAMPLE_N=inf;%6000; %%Number of cells to be sampled from each image
LOCCOLS={XCOORDLAB,YCOORDLAB};

sample_images=unique(double(slide_dm(:,'ImageNumber')))';  cnt=0;
fillIdx=0;
if ~exist('subI','var')
    subI=cell(1,length(sample_images));%%SUBSET of all cell's indices --> grouped by a cell's source image
    fillIdx=1;
end
% slide_dm=calculateDistanceToNeighbors(slide_dm,XCOORDLAB,YCOORDLAB);%%So that you can pick cells from dense regions
% figure('name','DTN')
% hist(slide_dm(:,'N_Neighbors'),50);
% return
imgI=cell(1,length(sample_images));%%ALL cell's indices --> grouped by a cell's source image
for i=sample_images
    cnt=cnt+1;
    i1=find(slide_dm(:,'ImageNumber')==i);% & slide_dm(:,'N_Neighbors')>=2);
    if fillIdx && length(i1)>SAMPLE_N*1.5
        idx1=unique(randi(length(i1),SAMPLE_N,1));%%sampling
        subI{cnt}=i1(idx1);
    elseif length(subI{cnt})>SAMPLE_N*1.5
        subI{cnt}=randsample(subI{cnt},SAMPLE_N);
    end
    imgI{cnt}=i1;
end

X=slide_dm(subI{1},LOCCOLS); %%img 1 - root image, to which other images must be adjusted - this one never changes coordinates!
% try
%     visualizeOverlayedImages(base_Image,slide_dm,imgI,XCOORDLAB,YCOORDLAB);%%Visualize before adjustments
% catch err
%     disp(err)
% end

%%Calculate adjustments for each image such that they match to first image
success=0;
adjustment=struct('rotationCenter_X',[],'rotationCenter_Y',[],'rotationAngle',[],'rotationDirection','','shift',[]);%%Overlay Adjustment parameters to be inferred
for imgIdx=2:length(subI)
    adjustment(imgIdx).rotationAngle=0;
    ITER_IMPROVEMENTS=nan(1,length(MAXITER));
    for iter=1:MAXITER
        Y=slide_dm(subI{imgIdx},LOCCOLS); %%reference position relative to which to shift
        %%%%%%%%%%%%%%%%%%%%%%%
        %%%%rotation center%%%%
        %         Dloc= pdist2(double(Y),double(X));
        %         diagIdx=-round(size(Y,1)/2):round(size(Y,1)/2);%% --> equivalent nuclei pairs expected along diagonal
        %         [~,ia]= min(nanmean(Dloc,2));%% --> location of those nuclei-pairs for which euclidean distance between figures is minimal
        %         xoyo=double(Y(ia,:)); %% rotation center
        xoyo=mean(Y);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%rotation angle THETA%%%%
        [~,~,THETA,adjustment(imgIdx)]=findAngleViaCorr(X,Y,xoyo, XCOORDLAB,YCOORDLAB,adjustment(imgIdx), iter);
        disp(['Rotation angle: ',num2str(THETA),' ',adjustment(imgIdx).rotationDirection])
        %         if iter>1 && THETA==0
        %             success=1;
        %             disp(['Overlay converged at iteration ',num2str(iter)])
        %             break;
        %         end
        
        %         Dloc= pdist2(double(Y_),double(Y));
        %         ITER_IMPROVEMENTS(iter)=diagSum(Dloc,diagIdx);
        %         Y=Y_;
        
        %%Save rotation and shift required to overlay this image with first
        adjustment(imgIdx).rotationCenter_X=xoyo(1);%rotation center
        adjustment(imgIdx).rotationCenter_Y=xoyo(2);%rotation center
        adjustment(imgIdx).rotationAngle=THETA;%rotation
        %%%%Apply inferred angle to overlay images%%%%
        [xr,yr] = rotateData(double(slide_dm(imgI{imgIdx},XCOORDLAB)),double(slide_dm(imgI{imgIdx},YCOORDLAB)),...
            adjustment(imgIdx).rotationCenter_X,adjustment(imgIdx).rotationCenter_Y,degtorad(adjustment(imgIdx).rotationAngle),adjustment(imgIdx).rotationDirection);
        slide_dm(imgI{imgIdx},LOCCOLS)=[xr;yr]';
        Y=slide_dm(subI{imgIdx},LOCCOLS);
        takeIfLast=slide_dm;
        
        %%Display overlay before this shift iteration
        figure('name', 'IPOX_HandE_Overlay','units','normalized','Position',[0,0,1,1],'Visible', 'on'); hold on;
        visi={'k.','r.','b.','c.'};
        for imgIdxI=1:2
            plot(double(slide_dm(imgI{imgIdxI},XCOORDLAB)),double(slide_dm(imgI{imgIdxI},YCOORDLAB)),visi{imgIdxI+1});
        end
        
        try
            xyAdjust=finetuneAlignment(X,Y,XCOORDLAB,YCOORDLAB,[0,1]);
        catch err
            disp(err)
            xyAdjust=[0,0];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%Shift direction using cross-correlation%%%%
        if iter>1
            allY=slide_dm(imgI{imgIdx},:);
            qs=quantile(slide_dm(slide_dm(:,'ImageNumber')==0,LOCCOLS),[0,1]);
            lastPercentDecline=inf;
            for prop=flip(0.7:0.05:0.9)-iter*1/100
                %         prop=0.9;
                innerX=slide_dm(slide_dm(:,'ImageNumber')==0,:);
                %             rX=unique(randi(size(innerX,1),3000000*prop,1));  rY=unique(randi(size(allY,1),3000000*prop,1));
                innerI=find(innerX(:,XCOORDLAB)>=mean(qs(:,1))-prop*mean(qs(:,1)) & innerX(:,XCOORDLAB)<=mean(qs(:,1))+prop*mean(qs(:,1)));
                innerI=intersect(innerI,find(innerX(:,YCOORDLAB)>=mean(qs(:,2))-prop*mean(qs(:,2)) & innerX(:,YCOORDLAB)<=mean(qs(:,2))+prop*mean(qs(:,2))));
                innerX=slide_dm(innerI,:);
                try
                    [xyAdjust_,percentDecline,dmHE_,dmIPOX_]=finetuneAlignment(innerX,allY,XCOORDLAB, YCOORDLAB,[0,1]);
                    if percentDecline<=55
                        slide_dm=[dmHE_;dmIPOX_]; % Perfect: done and done
                        success=1;
                        return;
                    elseif percentDecline<80 && percentDecline<lastPercentDecline
                        lastPercentDecline=percentDecline;
                        xyAdjust=xyAdjust_;
                        takeIfLast=[dmHE_;dmIPOX_];
                    end
                catch err
                    disp(err)
                end
            end
        end
        adjustment(imgIdx).xAdjust=xyAdjust(1);
        adjustment(imgIdx).yAdjust=xyAdjust(2);
        %         %%%%%%%%%%%%%%%%%%%%%%%
        %         %%%%Shift direction%%%%
        %         Dloc= pdist2(double(Y_),double(X)); Dloc=nanmin(Dloc); D=quantile(Dloc,0.85);
        %         shiftMat=1;         shiftMag=0;
        %         for theta=unique([0,(-1+iter):4:360])
        %             R = (50/iter)*repmat([1,1]*[cos(theta) -sin(theta); sin(theta) cos(theta)],size(Y_,1),1); %%direction
        %             y=Y_+R;
        %             Dloc= pdist2(double(y),double(X));
        %             Dloc=nanmin(Dloc);
        %             Dnew=quantile(Dloc,0.85);
        %             %             Dnew=diagSum(Dloc,diagIdx);
        %             if Dnew<D
        %                 D=Dnew;
        %                 shiftMat=R;
        %                 shiftMag=(2/iter);
        % %                 disp(['New best shift-direction: ',num2str(theta),' degrees']);
        %             end
        %         end
        
        %         %%%%%%%%%%%%%%%%%%%
        %         %%Shift magnitude%%
        %         D=inf;
        %         %         D= mean(quantile(Dloc,0.005));
        %         for s=unique([0,(-1+iter):5:300])
        %             y=Y_+shiftMat*s;
        %             Dloc= pdist2(double(y),double(X));
        %             Dloc=nanmin(Dloc);
        %             Dnew=quantile(Dloc,0.85);
        %             %             Dnew= mean(quantile(Dloc,0.005));
        % %             disp(['Shift-magnitude: ',num2str(s),'; DistanceM: ',num2str(Dnew)]);
        %             if Dnew<=D
        %                 D=Dnew;
        %                 shiftMag=s;
        %             elseif Dnew*0.95>=D
        %                 break;
        %             end
        %         end
        
        %         %%Combine direction & magnitude
        %         Y_=Y_+shiftMat*shiftMag;
        %         Y=Y_;
        
        
        %         adjustment(imgIdx).shift=shiftMat(1,:)*shiftMag;%shift
        if iter==MAXITER
            slide_dm=takeIfLast;
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%Apply inferred adjustment to overlay images%%%%
            x=slide_dm(imgI{imgIdx},LOCCOLS);
            %         x=x+repmat(adjustment(imgIdx).shift,size(x,1),1);
            x(:,1)=x(:,1)+adjustment(imgIdx).xAdjust;
            x(:,2)=x(:,2)+adjustment(imgIdx).yAdjust;
            slide_dm(imgI{imgIdx},LOCCOLS)=x;
            %%Quantify impact of adjustment on alignment
            qs=[0:0.05:1];
            ITER_IMPROVEMENTS(iter)=mean(mean(quantile(double(slide_dm(imgI{1},LOCCOLS)),qs)-quantile(x,qs)));
            disp(['Average alignment distance: ',num2str(ITER_IMPROVEMENTS(iter))]);
        end
        
    end
end

if ~success
    disp(['Overlay did not converge within ',num2str(MAXITER),' iterations.'])
end

% disp('Fine-tuning overlay using individual cell coordinates..')
% for imgIdx=2:length(subI)
%     [~, xyOffsets]=interSlideCoordinateCorrelation(slide_dm(subI{1},LOCCOLS), slide_dm(subI{imgIdx},LOCCOLS), XCOORDLAB, YCOORDLAB);
%     slide_dm(imgI{imgIdx},LOCCOLS)= slide_dm(imgI{imgIdx},LOCCOLS)+repmat(xyOffsets,length(imgI{imgIdx}),1);
% end

try
    visualizeOverlayedImages(base_Image,slide_dm,imgI,XCOORDLAB,YCOORDLAB);%%Visualize after adjustments
catch err
    disp(err)
end
end

function [dmIPOX,xyOffsets]=interSlideCoordinateCorrelation(dmHE, dmIPOX, xI, yI)
import bioma.data.*
%%Estimate inter-nuclei distances
tmp=randi(size(dmHE,1),100,1);
D=pdist2(double(dmHE(tmp,{xI,yI})), double(dmHE(setdiff(1:size(dmHE,1),tmp),{xI,yI})));
D=nanmin(D,[],2); D=nanmax(D);
NCELLS=4;
allC=[];
heIndices=cell(1,3);
for cnt=1:length(heIndices)
    %%NCELLS cells closest to center on H&E
    d=[0,0];
    while d(2)<=D %%pick cells that have a sparse neigborhood
        k=randi(size(dmHE,1),1,1);
        meanX=double(dmHE(k,xI));        meanY=double(dmHE(k,yI));
        x=pdist2([meanX,meanY],double(dmHE(:,{xI,yI})));
        [d,heI]=sort(x);
    end
    heI=heI(1:NCELLS);
    heIndices{cnt}=heI;
    
    %%NCELLS*2 cells closest to center on IPOX
    x=pdist2([meanX,meanY],double(dmIPOX(:,{xI,yI})));
    [~,ipoxI]=sort(x); ipoxI=ipoxI(1:ceil(NCELLS*2));
    
    %%Assess how well pattern of IPOX cells that matches H&E cell pattern
    ipoxPerm=combinator(length(ipoxI),length(heI),'p');
    c=DataMatrix(nan(size(ipoxPerm,1),8),cellstr(num2str(ipoxI(ipoxPerm))),{'rx','Px','ry','Py','r','P','VarX','VarY'});
    for i=1:size(ipoxPerm,1)
        if mod(i,200)==0
            disp(['Testing permutation ',num2str(i),' out of ',num2str(size(ipoxPerm,1)),' (iteration ',num2str(cnt),' out of ',num2str(length(heIndices)),')']);
        end
        perm=ipoxPerm(i,:);
        %%corr
        [rx,px]=corr(double(dmHE(heI,xI)),double(dmIPOX(ipoxI(perm),xI)));
        [ry,py]=corr(double(dmHE(heI,yI)),double(dmIPOX(ipoxI(perm),yI)));
        %%offset
        offSVar=var([dmHE(heI,xI)-dmIPOX(ipoxI(perm),xI), dmHE(heI,yI)-dmIPOX(ipoxI(perm),yI)]);
        
        c(i,:)=[rx,px,ry,py,mean([rx,ry]),mean([px,py]),offSVar];
    end
    c=addColumn(c,'Iter',cnt);
    allC=[allC;c];
end
[~,ia]=sort(mean(allC(:,{'VarX','VarY'}),2)); %[~,ia]=sort(double(allC(:,'P')));
allC=allC(ia,:);
%%Choose best match
ipoxI=str2num(char(regexp(char(rownames(allC(1,:))),'  ','split')'));
heI=heIndices{double(allC(1,'Iter'))};
xyOffsets=mean([dmHE(heI,xI)-dmIPOX(ipoxI,xI), dmHE(heI,yI)-dmIPOX(ipoxI,yI)]);
dmIPOX(:,{xI,yI})=dmIPOX(:,{xI,yI})+repmat(xyOffsets,size(dmIPOX,1),1);
end

function visualizeOverlayedImages(base_Image,slide_dm,imgI,XCOORDLAB,YCOORDLAB)

%%%%%%%%%%%%%%%
%%%Visualize%%%
% B=imread(base_Image);
% if size(B,3)>3
%     B=B(:,:,1:3);
% end
figure('name', 'IPOX_HandE_Overlay','units','normalized','Position',[0,0,1,1],'Visible', 'on'); hold on;
% imshow(B,'InitialMag','fit');
visi={'k.','r.','b.','c.'};
for imgIdx=1:length(imgI)
    plot(double(slide_dm(imgI{imgIdx},XCOORDLAB)),double(slide_dm(imgI{imgIdx}),YCOORDLAB),visi{imgIdx});
end

end


function dm=calculateDistanceToNeighbors(dm,xI,yI)
dm=addColumn(dm,'N_Neighbors',nan); dm(:,'N_Neighbors')=nan;
for i=1:size(dm,1)
    x=pdist2(double(dm(i,{xI,yI})),double(dm(:,{xI,yI})));
    nearI=setdiff(find(x<70),i);
    dm(i,'N_Neighbors')=length(nearI);
end
end

function [Y_,D,THETA,adjustment]=findAngleViaDistance(X,Yi,xoyo,XCOORDLAB,YCOORDLAB,adjustment)
import bioma.data.*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%rotation angle THETA%%%%
Y_=Yi;
adjustment.rotationDirection='clockwise'; %%Dummy default
THETA=0;       Dloc= pdist2(double(Yi),double(X));
D=mean(quantile(Dloc,0.01));%diagSum(Dloc,diagIdx);
for rotDirection={'clockwise','anticlockwise'}
    for theta=0:180
        [xr,yr] = rotateData(double(Yi(:,XCOORDLAB)),double(Yi(:,YCOORDLAB)),...
            xoyo(1),xoyo(2),degtorad(theta),char(rotDirection));
        y=DataMatrix([xr;yr]',rownames(Yi),{XCOORDLAB,YCOORDLAB}); %%Rotate
        Dloc= pdist2(double(y),double(X));
        Dnew=mean(quantile(Dloc,0.01));%diagSum(Dloc,diagIdx); %%choose angle that minimizes distance between equivalent nuclei pairs
        %         disp(['Angle: ',num2str(theta),'; Cummulative distance: ',num2str(Dnew)]);
        if Dnew<D
            D=Dnew;
            Y_=y;
            THETA=theta;
            adjustment.rotationDirection=char(rotDirection);
        elseif Dnew*0.8>=D
            break;
        end
    end
end
end

function [Y_,D,THETA,adjustment]=findAngleViaCorr(X,Yi,xoyo,XCOORDLAB,YCOORDLAB,adjustment,iter)
import bioma.data.*
if size(X,1)>7000
    X=X(unique(randi(size(X,1),7000,1)),:);
end
if size(Yi,1)>7000
    Yi=Yi(unique(randi(size(Yi,1),7000,1)),:);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%rotation angle THETA%%%%
% % % figure('name','CorrInSearchFor RightAngle'); hold on
% xO=struct('upper',[],'right',[]);%,'lower',[],'left',[]);
% for side=fieldnames(xO)'
%     xOutlineI=getSlideOutline(X, XCOORDLAB, YCOORDLAB,char(side),0.15);
%     yOutlineI=getSlideOutline(Yi, XCOORDLAB, YCOORDLAB,char(side),0.15);
%     tmpSz=min(length(yOutlineI),length(xOutlineI));
%     yI=unique(randi(length(yOutlineI),tmpSz*4,1),'stable');
%     xI=unique(randi(length(xOutlineI),tmpSz*4,1),'stable');
%     yI=yI(1:min(length(xI),length(yI))); xI=xI(1:min(length(xI),length(yI)));
% %         plot(X(xOutlineI(sort(xI)),XCOORDLAB), X(xOutlineI(sort(xI)),YCOORDLAB),'k*');
% %         plot(Yi(yOutlineI(sort(yI)),XCOORDLAB), Yi(yOutlineI(sort(yI)),YCOORDLAB),'r*');
%     xO=setfield(xO,char(side),[xOutlineI(sort(xI)); yOutlineI(sort(yI))]);
% end
% % % hold off;

rOutlineI=getSlideOutline(X, XCOORDLAB, YCOORDLAB,'right',0.15);
uOutlineI=getSlideOutline(X, XCOORDLAB, YCOORDLAB,'upper',0.15);
loOutlineI=getSlideOutline(X, XCOORDLAB, YCOORDLAB,'lower',0.15);
leOutlineI=getSlideOutline(X, XCOORDLAB, YCOORDLAB,'left',0.15);

Y_=Yi;
adjustment.rotationDirection='clockwise'; %%Dummy default
THETA=0;      D=[inf,inf];
%         D=mean(quantile(Dloc,0.01));%diagSum(Dloc,diagIdx);
for rotDirection={'clockwise','anticlockwise'}
    for theta=unique([0,(-0.1+0.1*iter):(2/iter):180])
        [xr,yr] = rotateData(double(Yi(:,XCOORDLAB)),double(Yi(:,YCOORDLAB)),...
            xoyo(1),xoyo(2),degtorad(theta),char(rotDirection));
        y=DataMatrix([xr;yr]',rownames(Yi),{XCOORDLAB,YCOORDLAB}); %%Rotate
        yuOutlineI=getSlideOutline(y, XCOORDLAB, YCOORDLAB,'upper',0.15);
        
        %%Alternative 1
        %         rX=nan(size(X,1),1);
        %         for i=[xOutlineI,yOutlineI]
        %             diff=sum(bsxfun(@minus, double(y(:,{XCOORDLAB,YCOORDLAB})), double(X(i,{XCOORDLAB,YCOORDLAB}))),2);
        %             [~,ia]=nanmin(abs(diff));
        %             rX(i)=diff(ia);
        %         end
        %%Alternative 2
        %         rX=pdist2(double(X(:,{XCOORDLAB,YCOORDLAB})),double(y(:,{XCOORDLAB,YCOORDLAB})));
        %         rX=min(rX);
        %         Dnew=std(rX(uOutlineI))*length(uOutlineI)+std(rX(rOutlineI))*length(rOutlineI)+...
        %             std(rX(loOutlineI))*length(loOutlineI)+std(rX(leOutlineI))*length(leOutlineI);
        %%Alternative 3 (ORIG)
        rX=pdist2(double(X(:,{XCOORDLAB,YCOORDLAB})),double(y(:,{XCOORDLAB,YCOORDLAB})));
        rX=min(rX);
        Dnew1=mean(rX);
%         Dnew2=Dnew1;
        %%Alternative 4
        crr=xcorr2((double(X(uOutlineI,{XCOORDLAB,YCOORDLAB}))), (double(y(yuOutlineI,{XCOORDLAB,YCOORDLAB}))));
        [~,snd] = max(crr(:)); crrsub=crr(round(snd*0.85):round(snd*1.15));
        Dnew2=-100*max(crrsub([1:1000,(length(crrsub)-1000):length(crrsub)]))/max(crrsub);
      
        if all([Dnew1,Dnew2]<D)
            D=[Dnew1,Dnew2];
            Y_=y;
            THETA=theta;
            adjustment.rotationDirection=char(rotDirection);
        end
    end
end

% figure('name','CorrInSearchFor RightAngle'); hold on
% plot(X(xO.upper(1,:),XCOORDLAB), X(xO.upper(1,:),YCOORDLAB),'k*');
% plot(X(xO.right(1,:),XCOORDLAB), X(xO.right(1,:),YCOORDLAB),'k*');
% plot(Y_(xO.upper(2,:),XCOORDLAB), Y_(xO.upper(2,:),YCOORDLAB),'r*');
% plot(Y_(xO.right(2,:),XCOORDLAB), Y_(xO.right(2,:),YCOORDLAB),'r*');
% hold off
end

%%%%rotation center%%%%
% figure('name','RotationCenterChoice')
% subplot(1,2,1)
% plot(X(:,LOCCOLS{1}),nanmean(Dloc,2),'r*')
% xlabel('X-coordinate');ylabel('mean pairwise distance')
%
% subplot(1,2,2)
% plot(nanmean(Dloc,2),X(:,LOCCOLS{2}),'r*')
% ylabel('Y-coordinate');xlabel('mean pairwise distance')
%
