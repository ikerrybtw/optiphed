function net=trainClassifier(sampleName,path2features)

import bioma.data.*

global MARKERF
global XCOORDLAB
global YCOORDLAB
global COLORRES
XCOORDLAB='AreaShape_Center_X';
YCOORDLAB='AreaShape_Center_Y';
COLORRES=200;
PATH2HYPERPARAMETERS='data/NeuralNetworkHyperparameters.xlsx';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Read names of features used for classification%%%%
cassFeatures=importdata(path2features);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Read values of features calculated from image for each nucleus%%%%
image1TIF=['data/training/',sampleName,'.tif'];
hande_dm=readNucleiFeatures(sampleName,'_AllCols.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Display training region%%%%
figure('name', char(sampleName),'units','normalized','Position',[0,0,0.75,1],'Visible','on');
try
    B=imread(image1TIF);
    if  length(size(B))>2
        B=B(:,:,1:3);
    end
    subplot(2,3,1); image(B);%,'InitialMag','fit');
catch
    disp(['No image found for ',char(sampleName)])
end
subplot(2,3,2); 
drawNucleiOutlines(hande_dm,MARKERF,COLORRES,2); %%[XCOORDLAB,'_Orig'], [YCOORDLAB,'_Orig']
set(gca,'Ydir','reverse')
set(gca,'Color',[0.8 0.8 0.8]);

%%Group into IPOX+(tumor cells) vs. IPOX-(normal cells)
subplot(2,3,4); hist(double(hande_dm(:,MARKERF)),100);
cl=nan(size(hande_dm,1),1); cl(hande_dm(:,MARKERF)>0.92)=1; cl(hande_dm(:,MARKERF)<0.75)=0;
cl=imputeClassFromNeigbors(hande_dm,0.005,XCOORDLAB,YCOORDLAB,cl);
subplot(2,3,3); plot(double(hande_dm(cl==1,[XCOORDLAB,'_Orig'])), double(hande_dm(cl==1,[YCOORDLAB,'_Orig'])),'k*'); hold on
subplot(2,3,3); plot(double(hande_dm(cl==0,[XCOORDLAB,'_Orig'])), double(hande_dm(cl==0,[YCOORDLAB,'_Orig'])),'cyan*'); set(gca,'Ydir','reverse')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Train neural network classifier%%%%
X=double(hande_dm(~isnan(cl),cassFeatures));      t=cl(~isnan(cl));
fI=find(all(isfinite(X),2)); X=X(fI,:);  t=t(fI);
disp(['Using ',num2str(size(X,2)),' features to train classifier on ',num2str(length(fI)),' nuclei'])
[net, hyperparams] = chooseNetwork( X, t,  PATH2HYPERPARAMETERS);
disp(hyperparams)
    

end