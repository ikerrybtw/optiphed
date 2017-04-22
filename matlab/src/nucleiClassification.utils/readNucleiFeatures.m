function hande_dm=readNucleiFeatures(sampleName,suffix)
import bioma.data.*

image1F=['data/training/',sampleName,suffix];
A = importdata(image1F);
hande_dm=DataMatrix(A.data(:,2:size(A.data,2)),1:size(A.data,1),strsplit(char(strtrim(cellstr(A.textdata))),'\t'));
disp(['Read ',num2str(size(hande_dm,2)),' features for ',num2str(size(hande_dm,1)),' nuclei'])
%%Exclude spurious nuclei with very small area --> they are probably artifacts
if any(strcmp(colnames(hande_dm),'AreaShape_Area'))
    ii=find(hande_dm(:,'AreaShape_Area')<=2*min(hande_dm(:,'AreaShape_Area')));
    disp(['Removing ',num2str(length(ii)),' spurious nuclei (',num2str(round(100*length(ii)/size(hande_dm,1),2)),'%)'])
    hande_dm(ii,:)=[];
end
end