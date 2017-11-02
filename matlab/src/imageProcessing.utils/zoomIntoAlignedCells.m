function regions=zoomIntoAlignedCells(fig, dm, cols, clusterCols ,K)
% K=round(size(dm,1)/NCELLS);
disp(['Looking for ',num2str(K),' clusters...'])
Z=linkage(double(dm(:,clusterCols)),'average','euclidean');
cl = cluster(Z,'maxclust',K);     
% cl=kmeans(double(dm(:,clusterCols)),K);
clMeans=grpstats(double(dm(:,clusterCols)),cl,'mean');
cl=kmeans(double(dm(:,clusterCols)),K,'start',sort(clMeans)); %%Ensure that clusters are in same order
figure(fig)
oldFigName=get(fig,'name');
xL=get(gca,'XLim');  yL=get(gca,'YLim');

regions={};
for c= unique(cl)' 
    ii=find(cl==c);
    xyL=quantile(dm(ii,cols),[0,1]);
    xyL(1,:)=xyL(1,:)*0.99;  
    xyL(2,:)=xyL(2,:)*1.01;  
    set(gca,'XLim',xyL(:,1)','YLim',xyL(:,2)')
    export_fig([oldFigName,'_Cluster',num2str(c)],'-tif');%,'-nocrop')
    %%Save mat files too
    regions=[regions,[oldFigName,'_Cluster',num2str(c),'.txt']];
    dmwrite(dm(ii,:),regions{length(regions)})
end
end