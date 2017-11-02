function overlayIPOXonHandE(indir,sampleName)
%% Directory 'indir' must contain two subdirectories: 1. HandE; 2. IPOX. Each subdirectory must contain a folder starting with the sample ID
%% Global and local alignment of H&E and IPOX stained nuclei identified by CellProfiler

%indir='/mnt/ix2/nandor/Projects/ITHinferenceFromHandE/data/PhenotypicSPsFromHandE_trainingSet/IndividualNuclei_TrainingSet/Nuclei/UCLA_round2/layer2'
%samples={'1957_T_9668','2019_T_9717','2011_T_9699','1939_T_9661','1932_T_9876', '1928_T_9644','1792_T_9500', '1915_T_10464', '1910_T_9637', '1794_T_9504','1919_T_10779'};
% for sample = flip(samples)
% disp(['Processing ',char(sample),'...'])
% overlayIPOXonHandE(indir,char(sample));
% close all;
%  end

import bioma.data.*
%root2='/Users/noemi/Desktop/';       root='/Users/noemi/';
root='/mnt/ix2/nandor/'; root2=root;
prjRoot=[root,'Projects/ITHinferenceFromHandE/'];
addpath([root2,'/Functionalities/Skripts/'])
addpath([root2,'/Functionalities/Skripts/export_fig/'])
addpath([root,'/Projects/ITHinferenceFromHandE/code/imageProcessing/'])
addpath([root,'/Projects/code/PMO_matlab/src_matlab/BioSample_PhenotypicSubpopulation/'])

sampleID=strsplit(sampleName,'_');
sampleID=char(strcat(strjoin(strcat(sampleID(1:2),'_'),''),sampleID(3)));

NCELLS=500;
COLORRES=200;
MARKERF='Intensity_MaxIntensityEdge_Costum'; %%The IPOX staining marker column
MARKERS={'HandE','IPOX'};
XCOORDLAB='AreaShape_Center_X';
YCOORDLAB='AreaShape_Center_Y';
IPOXONLYCOLS={MARKERF,'Intensity_UpperQuartileIntensity_Costum'};
QCSTATUSFILE=struct('global','QUALITYCHECK_GLOBALALIGNMENT.txt','local','QUALITYCHECK_LOCALALIGNMENT.txt'); %%to be changed by user after quality is manually confirmed OK
[~,tmp]=fileparts([indir,filesep,MARKERS{1},filesep,sampleName]);
LAYER=strsplit(tmp,'layer');
LAYER=['layer',LAYER{2}];

%%Outdir structure per tiers
tier1Out=strrep(indir,'/Nuclei/','/Overlay/Tier1/');
tier2Out=strrep(indir,'/Nuclei/','/Overlay/Tier2/');
tier3Out=strrep(indir,'/Nuclei/','/Overlay/Tier3/');
tier4Out=strrep(indir,'/Nuclei/','/Overlay/Tier4/');
for outD={tier1Out,tier2Out,tier3Out,tier4Out}
    if ~exist(char(outD),'dir')
        mkdir(char(outD))
    end
end
tier1MatBeforeOverlay=[tier1Out,filesep,sampleID,'_',LAYER,'.mat'];
tier1Mat=[tier1Out,filesep,sampleID,'_',LAYER,'_overlayedAdjacentSlides.mat'];
tier3Mat=[tier3Out,filesep,sampleID,'_',LAYER,'_IPOX_to_HandE_AssignedNuclei.mat'];

%%Informs user about new file
info_=[];
if exist(QCSTATUSFILE.global,'file')
    info_=DataMatrix('File',QCSTATUSFILE.global); info_=DataMatrix(info_,rownames(info_),strtrim(colnames(info_)));
end
overlayQNotOK=(~isempty(info_) && any(strcmp(rownames(info_),sampleID)) && info_(sampleID,'qualityOK')<0);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%READING NUCLEI FILES%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Reading H&E and IPOX Nuclei files for ', sampleID]);
dm=[]; tifs=struct(); loaded=0;
if exist(tier1MatBeforeOverlay,'file')
    disp(['Loading previously saved IPOX and H&E matrix for ', sampleID,' from ',tier1MatBeforeOverlay]);
    load(tier1MatBeforeOverlay); %%overloads dm
    loaded=1;
end
for i=1:length(MARKERS)
    tmp=dir([indir,filesep,MARKERS{i},filesep,sampleID,'*',LAYER]);
    iD=find([tmp.isdir]);
    nucleiFile=[indir,filesep,MARKERS{i},filesep,tmp(iD).name,filesep,'cellProfilerOut',filesep,'Nuclei.txt'];
    tifs=setfield(tifs,MARKERS{i},[indir,filesep,MARKERS{i},filesep,tmp(iD).name,filesep,tmp(iD).name,'.tif']);
    if ~loaded
        if (~exist(tier1Mat,'file') && ~exist(tier3Mat,'file') ) || overlayQNotOK
            disp(['Reading ', nucleiFile,'...']);
            
            addCols={};
            if i>1
                addCols=IPOXONLYCOLS;
            end
            y=DataMatrix('File',nucleiFile,'Columns', [XCOORDLAB,YCOORDLAB,addCols]);
            y=addColumn(y,'ImageNumber',i-1);
            
            dm= joinMatrices( dm,y);
            save(tier1MatBeforeOverlay,'dm');
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%GLOBAL ALIGNMENT%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (~exist(tier1Mat,'file') && ~exist(tier3Mat,'file') ) || overlayQNotOK
    disp(['Performing global alignment for ', sampleID]);
    %%Save original XY-coordinates
    dm=addColumn(dm,[XCOORDLAB,'_Orig'],nan);  dm=addColumn(dm,[YCOORDLAB,'_Orig'],nan);
    dm(:,{[XCOORDLAB,'_Orig'],[YCOORDLAB,'_Orig']})=dm(:,{XCOORDLAB,YCOORDLAB});
    %%Remove small chuncks of tissue
    firstR=1; toRm=[];
    while firstR || ~isempty(toRm)
        toRm=[]; firstR=0;
        dm=removeNoiseFromImage(dm,XCOORDLAB,YCOORDLAB);%%Remove noise
        for img=unique(double(dm(:,'ImageNumber')))'
            for coordAndSign={{XCOORDLAB,'-'},{XCOORDLAB,'+'},{YCOORDLAB,'-'},{YCOORDLAB,'+'}}
                cS=coordAndSign{1}; COORD=cS{1};
                if strcmp(cS{2},'-')
                    dm(:,COORD)=-1*dm(:,COORD); %%Reverse sign to look at opposite end
                end
                rangeI=max(dm(:,COORD)) - min(dm(:,COORD));
                dmI=find(dm(:,'ImageNumber')==img & dm(:,COORD)<min(dm(:,COORD))+0.3*rangeI);
                obj=gmdistribution.fit(double(dm(dmI,{XCOORDLAB,YCOORDLAB})),2);
                cl=obj.cluster(double(dm(dmI,{XCOORDLAB,YCOORDLAB})));
%                 [~,clI]=min(obj.ComponentProportion); %%The smaller chunck 
                
                [medX,stdX]=grpstats(double(dm(dmI,COORD)),cl,{'median','std'});
                [~,clI]=min(medX); %%The leftmost/lowermost chunck
                [~,clJ]=max(medX); %%The rightmost/uppermost chunck
                stDev=std(dm(dmI(cl==clI),{XCOORDLAB,YCOORDLAB}));
                rI=unique(randi(length(dmI),1000,1)); %%Sample index to calc distance btw. clusters
                d=pdist2(double(dm(dmI(intersect(rI,find(cl==clI))),{XCOORDLAB,YCOORDLAB}))   ,    double(dm(dmI(intersect(rI,find(cl==clJ))),{XCOORDLAB,YCOORDLAB})));
                clDist=quantile(reshape(d,prod(size(d)),1),0.0); %%Cluster distance metric
                if clDist>200 && pdist(obj.mu)>pdist([obj.mu(clI,:);obj.mu(clI,:)+stDev]) %%If distance between peaks > distance between 1st peak and its std
                    toRm=[toRm; dmI(cl==clI)];
                end
                dm(:,COORD)=abs(dm(:,COORD));
            end
        end
        if ~isempty(toRm)
            disp(['Removing small tissue chunck consisting of ',num2str(length(toRm)),' cells']);
            dm(toRm,:)=[];
        end
    end
    close all; plotOverlay(dm,XCOORDLAB,YCOORDLAB,tier1Out,sampleID,LAYER);
    [dmO,keepI0,keepI1, dropI0, dropI1]=getSlideOutlines(dm, XCOORDLAB, YCOORDLAB);%%Just the image outlines
    dm=overlayAdjacentSlides(dm,'', XCOORDLAB, YCOORDLAB,{keepI0,keepI1}, 5);
    for fig={'HandE_On_IPOXbackground_IPOXQ_0-1.fig','CrossCorrMax_0-1.fig'}
        if exist(char(fig),'file')
            movefile(char(fig), [tier1Out,filesep,sampleID,'_',char(fig)]);
        end
    end
    % %%Split figure in 4 quartiles when overlaying
    % dm_=dm;
    % tmp=quantile(dm(:,XCOORDLAB),[0,1]); xL=tmp(1):(tmp(2)-tmp(1))/2:tmp(2);
    % tmp=quantile(dm(:,YCOORDLAB),[0,1]); yL=tmp(1):(tmp(2)-tmp(1))/2:tmp(2);
    % for qx=1:(length(xL)-1)
    %     ii=find(dm(:,XCOORDLAB)>=xL(qx) & dm(:,XCOORDLAB)<=xL(qx+1));
    %     for qy=1:(length(yL)-1)
    %         qName=['Qx',num2str(floor(qx)),'Qy',num2str(floor(qy))];
    %         ij=find(dm(:,YCOORDLAB)>=yL(qy) & dm(:,YCOORDLAB)<=yL(qy+1));
    %         qI=intersect(ii,ij);
    %         disp([num2str(length(qI)),' nuclei in quartile ',qName])
    %         [~,i0]=intersect(qI,keepI0);        [~,i1]=intersect(qI,keepI1);
    %         tmp=overlayAdjacentSlides(dm(qI,:), [], XCOORDLAB, YCOORDLAB, {i0,i1}, 4);
    %         dm_(qI,:)=tmp;
    %         set(gcf,'name',[sampleID,'_',LAYER,'_',get(gcf,'name'),'_AllIntensities_',qName]); title(sampleID)
    %     end
    % end
    % dm=dm_;
    close all
    save(tier1Mat,'dm');
    %%Assign closest nuclei even at this lower overlay resolution as intermediate result
    A=dm(dm(:,'ImageNumber')==0,:);    B=dm(dm(:,'ImageNumber')>0,:); 
    A=assignClosestNeighbor(A,B,XCOORDLAB,YCOORDLAB, {'Intensity_MaxIntensityEdge_Costum','Intensity_UpperQuartileIntensity_Costum'}); 
    dmwrite(A,strrep(tier1Mat,'overlayedAdjacentSlides.mat',[MARKERS{1},'_ClusterAll.txt']) )

    %%Informs user about new file
    if exist(QCSTATUSFILE.global,'file')
        info_=DataMatrix('File',QCSTATUSFILE.global); info_=DataMatrix(info_,rownames(info_),strtrim(colnames(info_)));
    end
    info=DataMatrix(zeros(1,1),{sampleID},{'qualityOK'});
    info=[info_;info];
    [~,ia]=unique(rownames(info)); info=info(ia,:);
    info(sampleID,'qualityOK')=0;
    dmwrite(info,QCSTATUSFILE.global, 'APPEND',0);
elseif exist(tier1Mat,'file')
    disp(['Loading global alignment for ', sampleID,' from ',tier1Mat]);
    load(tier1Mat)
end
%%(re-)plot global overlay
close all
figure('name', 'IPOX_HandE_Overlay','units','normalized','Position',[0,0,3,2],'Visible', 'off');
visi={'k.','r.','b.','c.'};
for img=unique(double(dm(:,'ImageNumber')))'
    plot(double(dm(dm(:,'ImageNumber')==img,XCOORDLAB)),double(dm(dm(:,'ImageNumber')==img,YCOORDLAB)),visi{img+1}); hold on
end
set(gcf,'name',[tier1Out,filesep,sampleID,'_',LAYER,'_',get(gcf,'name'),'_AllIntensities']); title(sampleID);
set(gca,'XLim',[min(dm(:,XCOORDLAB)),max(dm(:,XCOORDLAB))],'YLim',[min(dm(:,YCOORDLAB)),max(dm(:,YCOORDLAB))])
savefigs();
%%
info=DataMatrix('File',QCSTATUSFILE.global); info=DataMatrix(info,rownames(info),strtrim(colnames(info)));

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%LOCAL ALIGNMENT%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Only performed once QC of global alignment is confirmed OK by user
while info(sampleID,'qualityOK')==0
    info=DataMatrix('File',QCSTATUSFILE.global); info=DataMatrix(info,rownames(info),strtrim(colnames(info)));
    disp(['******Check global overlay quality under: ',get(gcf,'name')])
    disp(['If QC is OK, set ', sampleID, ' in ',QCSTATUSFILE.global,' to 1 to proceed with local alignment, otherwise set to -1 if QC failed'])
    pause(20);
end
if info(sampleID,'qualityOK')==-1
    return;
end

%%Informs user about new file
localOverlayQNotOK=boolean(0);
if exist(QCSTATUSFILE.local,'file')
    info_=DataMatrix('File',QCSTATUSFILE.local); info_=DataMatrix(info_,rownames(info_),strtrim(colnames(info_)));
    ii=find(cellfun(@(x) ~isempty(x),strfind(rownames(info_),sampleID)));
    localOverlayQNotOK = ~isempty(ii) && all(info_(ii,'qualityOK')<0);
end

if ~exist(tier3Mat,'file') || overlayQNotOK || localOverlayQNotOK
    %%Progress to next tier
    copyfile([tier1Out,filesep,sampleID,'*'],[tier2Out,filesep]);
    %%Save global alignment coordinates before local adjustment
    dm=addColumn(dm,[YCOORDLAB,'_NoLocalAdjustment'],nan); dm=addColumn(dm,[XCOORDLAB,'_NoLocalAdjustment'],nan);
    dm(:,{[XCOORDLAB,'_NoLocalAdjustment'],[YCOORDLAB,'_NoLocalAdjustment']})=dm(:,{XCOORDLAB,YCOORDLAB});
    hande_dm=dm(dm(:,'ImageNumber')==0,:);          ipox_dm=dm(dm(:,'ImageNumber')~=0,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% Normalized coordinates%%%%
    LOCCOLSTONORM={XCOORDLAB,YCOORDLAB,[XCOORDLAB,'_NoLocalAdjustment'],[YCOORDLAB,'_NoLocalAdjustment']};
    ipox_dm(:,LOCCOLSTONORM)=normalize01(double(ipox_dm(:,LOCCOLSTONORM)));
    hande_dm(:,LOCCOLSTONORM)=normalize01(double(hande_dm(:,LOCCOLSTONORM)));
    
    %%%Assign IPOX nuclei to H&E nuclei%%%%
    disp(['Performing local alignment for ', sampleID]);
    done=[]; heI=[]; ipoxI=[];
    qs=setdiff(unique(0.5+randi(5E9,1,350)/(5E9*2),'stable'),done,'stable') ;
    for   q=qs
        [success,hande_dm_,ipox_dm_]=calculateDistanceOfHandENucleiToIPOX(hande_dm, ipox_dm,XCOORDLAB,YCOORDLAB, NCELLS,q);
        close all;
        if success
            %%Get index in original matrix
            [~,~,heI_]=intersect(double(hande_dm_(:,{[XCOORDLAB,'_NoLocalAdjustment'],[YCOORDLAB,'_NoLocalAdjustment']})), double(hande_dm(:,{[XCOORDLAB,'_NoLocalAdjustment'],[YCOORDLAB,'_NoLocalAdjustment']})),'rows','stable');
            [~,~,ipoxI_]=intersect(double(ipox_dm_(:,{[XCOORDLAB,'_NoLocalAdjustment'],[YCOORDLAB,'_NoLocalAdjustment']})), double(ipox_dm(:,{[XCOORDLAB,'_NoLocalAdjustment'],[YCOORDLAB,'_NoLocalAdjustment']})),'rows','stable');
            
            IPOXONLYCOLS=union(IPOXONLYCOLS, setdiff(colnames(hande_dm_),colnames(hande_dm)));
            %%add columns if they don't exist
            for c =setdiff(IPOXONLYCOLS,colnames(hande_dm))
                hande_dm=addColumn(hande_dm,char(c),nan);
            end
            %%Set values for the matching indices
            countMarked=sum(~isnan(hande_dm(:,MARKERF)));
            for c =IPOXONLYCOLS
                hande_dm(heI_,c)=hande_dm_(:,c);
            end
            countMarkedNew=sum(~isnan(hande_dm(:,MARKERF)));
            heI=unique([heI;heI_]);     ipoxI=unique([ipoxI;ipoxI_]);
            %%Add as new cluster only if sufficiently new values have been
            %%added
            disp([MARKERF,' info added for ',num2str(countMarkedNew-countMarked),' additional H&E nuclei'])
            if countMarkedNew-countMarked>=0.99*size(hande_dm_,1) %%99% of the H&E nuclei-info is new
                done=[done,q];
            end
        end
        if length(heI)>=150000
            break;
        end
    end
    assigned_hande=struct('dmHE',hande_dm,'heI',sort(heI),'dmIPOX',ipox_dm,'ipoxI',sort(ipoxI),'done',done);
    save(tier3Mat,'assigned_hande');
else
    disp(['Loading local alignment for ', sampleID,' from ',tier3Mat]);
    load(tier3Mat);
end

%%Visualize locally aligned regions
disp(['Visualizing locally aligned regions for ', sampleID]);
ipox_dm=assigned_hande.dmIPOX(assigned_hande.ipoxI,:);
hande_dm=assigned_hande.dmHE(assigned_hande.heI,:);
hande_dm(hande_dm(:,'ClosestIPOXNuclei')>exp(2),MARKERF)=0; %%Nuclei that could not be identified in IPOX as well are negative for staining.
%draw H&E
figHE=overlayOutlinesOnHandE(fileparts(tifs.HandE),strrep(getfield(dir(tifs.HandE),'name'),'.tif',''),'on',[0,0,3,2]);
hold on; drawNucleiOutlines(hande_dm,[XCOORDLAB,'_Orig'], [YCOORDLAB,'_Orig'],'AreaShape_MajorAxisLength','AreaShape_MinorAxisLength', 'AreaShape_Orientation', MARKERF,COLORRES,5);
set(figHE,'name',[tier3Out,filesep,sampleID,'_',LAYER,'_HandE_nucleiOverlayed']); title(sampleID);
title(sampleID);
savefigs();
regions=zoomIntoAlignedCells(figHE, hande_dm,{[XCOORDLAB,'_Orig'],[YCOORDLAB,'_Orig']},{XCOORDLAB,YCOORDLAB},  length(assigned_hande.done));%Zoom in & save locally aligned regions
regions=strrep(regions,'_HandE_nucleiOverlayed','_*_nucleiOverlayed'); %%To be saved later on: Include IPOX as well
close all
%%draw IPOX in same region
figIPOX=overlayOutlinesOnHandE(fileparts(tifs.IPOX),strrep(getfield(dir(tifs.IPOX),'name'),'.tif',''),'on',[0,0,3,2]);
% tmp=pdist2(double(ipox_dm(:,{XCOORDLAB,YCOORDLAB})) ,  double(hande_dm(:,{XCOORDLAB,YCOORDLAB})) );
% [~,ia]=min(tmp);
hold on; drawNucleiOutlines(ipox_dm,[XCOORDLAB,'_Orig'], [YCOORDLAB,'_Orig'],'AreaShape_MajorAxisLength','AreaShape_MinorAxisLength', 'AreaShape_Orientation', MARKERF,COLORRES,5);
set(figIPOX,'name',[tier3Out,filesep,sampleID,'_',LAYER,'_IPOX_nucleiOverlayed']); title(sampleID);
savefigs();
zoomIntoAlignedCells(figIPOX, ipox_dm,{[XCOORDLAB,'_Orig'],[YCOORDLAB,'_Orig']},{XCOORDLAB,YCOORDLAB}, length(assigned_hande.done));%Zoom in & save locally aligned regions
close all
%%Informs user about new file

if exist(QCSTATUSFILE.local,'file')
    info=DataMatrix('File',QCSTATUSFILE.local); info=DataMatrix(info,rownames(info),strtrim(colnames(info)));
    info=[info;DataMatrix(zeros(length(regions),1),regions,{'qualityOK'})];
    [~,ia]=unique(rownames(info)); info=info(ia,:);
    info(regions,'qualityOK')=0;
else
    info=DataMatrix(zeros(length(regions),1),regions,{'qualityOK'});
end
dmwrite(info,QCSTATUSFILE.local, 'APPEND',0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Move succesfull local re-alignment into last tier (4)%%%%
while any(info(regions,'qualityOK')==0)
    info=DataMatrix('File',QCSTATUSFILE.local); info=DataMatrix(info,rownames(info),strtrim(colnames(info)));
    disp(['******Check local overlay quality for each region under: ',tier3Out,filesep,sampleID,'*_nucleiOverlayed_Cluster*.tif'])
    disp(['If QC is OK, set corresponding region in ',QCSTATUSFILE.local,' to 1 to proceed to tier 4, otherwise set to -1 if QC failed'])
    for region=regions
        if info(region,'qualityOK')>=1 && length(dir(char(region)))
            disp(['Moving ', char(region),' to ' ,tier4Out,filesep])
            movefile(char(strrep(region,'.txt','.tif')),[tier4Out,filesep]);%%Progress to next tier
            movefile(char(region),[tier4Out,filesep]);%%Progress to next tier
        end
    end
    pause(20);
end
dmwrite(info(info(:,'qualityOK')>=1,:),[tier4Out,filesep,'LOCALALIGNMENT_TrainingRegions.txt'])

end




function plotOverlay(dm,XCOORDLAB,YCOORDLAB,tier1Out,sampleID,LAYER)
figure('name', 'IPOX_HandE_Overlay','units','normalized','Position',[0,0,3,2],'Visible', 'off');
visi={'k.','r.','b.','c.'};
for img=unique(double(dm(:,'ImageNumber')))'
    plot(double(dm(dm(:,'ImageNumber')==img,XCOORDLAB)),double(dm(dm(:,'ImageNumber')==img,YCOORDLAB)),visi{img+1}); hold on
end
set(gcf,'name',[tier1Out,filesep,sampleID,'_',LAYER,'_',get(gcf,'name'),'_AllIntensities']); title(sampleID);
set(gca,'XLim',[min(dm(:,XCOORDLAB)),max(dm(:,XCOORDLAB))],'YLim',[min(dm(:,YCOORDLAB)),max(dm(:,YCOORDLAB))])
savefigs();
end


