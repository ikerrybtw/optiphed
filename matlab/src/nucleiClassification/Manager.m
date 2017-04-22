addpath matlab/src/nucleiClassification
addpath matlab/src/nucleiClassification.utils
addpath matlab/src/utils/
PATH2FEATURES='data/cassette1Features.txt';
SAMPLE1='1794_T_9504_layer2_HandE_nucleiOverlayed_Cluster1';
SAMPLE2='2038_T_11986_layer2_HandE_nucleiOverlayed_Cluster5';
global MARKERF
MARKERF='Intensity_MaxIntensityEdge_Costum'; %%The IPOX staining marker column name, holding the gold standart classification

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%train classifier on image 1%%
net=trainClassifier(SAMPLE1,PATH2FEATURES);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%apply classifier on image1 and image2%%
cassFeatures=importdata(PATH2FEATURES);
for sample={SAMPLE1,SAMPLE2}
    dm=readNucleiFeatures(char(sample),'_AllCols.txt');
    dm=addColumn(dm,'Class',NaN);
    dm(:,'Class')=net(double(dm(:,cassFeatures))')';
    figure('name',[char(sample),'_Class'],'units','normalized','Position',[0,0,0.75,1],'Visible','on');
    drawNucleiOutlines(dm,'Class',200,2); set(gca,'Ydir','reverse')
    %%Assess performance:
    [s,v,th,au]=perfcurve(double(dm(:,MARKERF)>0.92), double(dm(:,'Class')>0.9),'1','XCrit','sens','YCrit','spec');
    disp([char(sample),': AUC=',num2str(au)])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%compare to gold standard%%
for sample={SAMPLE1,SAMPLE2}
    dm=readNucleiFeatures(strrep(char(sample),'HandE','IPOX'),'.txt');
    figure('name',[char(sample),'_IPOX'],'units','normalized','Position',[0,0,0.75,1],'Visible','on');
    drawNucleiOutlines(dm,MARKERF,200,2); set(gca,'Ydir','reverse');
end
