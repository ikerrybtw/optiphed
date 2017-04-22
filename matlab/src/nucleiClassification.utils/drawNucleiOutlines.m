function drawNucleiOutlines(hande_dm,intensityCol,colorResolution,msize)
global XCOORDLAB
global YCOORDLAB

if ~exist('msize','var')
    msize=1;
end

iK=find(~isnan(hande_dm(:,intensityCol)));
c=double(ceil(hande_dm(iK,intensityCol)*colorResolution));
cut=round(nanmin(hande_dm(:,intensityCol)*colorResolution))-1;
colors=flip(hot(1+round(nanmax(hande_dm(:,intensityCol)*colorResolution))-cut));


%%Now plot cells
hold on;
scatter(double(hande_dm(iK,XCOORDLAB)),double(hande_dm(iK,YCOORDLAB)),msize,colors(c-cut,:),'LineWidth',2);

xL=quantile(hande_dm(:,XCOORDLAB),[0,1]);
yL=quantile(hande_dm(:,YCOORDLAB),[0,1]);
set(gca,'XLim',xL','YLim',yL');
colorbar;
colormap(colors);

end