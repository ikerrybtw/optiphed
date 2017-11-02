function [dmO,keepI0,keepI1, dropI0, dropI1]=getSlideOutlines(dm, xI, yI)
[dm0,dropI0]=getSlideOutlinesPerImage(dm(dm(:,'ImageNumber')==0,:),xI,yI);
[dm1,dropI1]=getSlideOutlinesPerImage(dm(dm(:,'ImageNumber')~=0,:),xI,yI);
dmO=[dm0;dm1];
[~,keepI0]=intersect(double(dm(:,{xI,yI,'ImageNumber'})),double(dm0(:,{xI,yI,'ImageNumber'})),'rows');
[~,keepI1]=intersect(double(dm(:,{xI,yI,'ImageNumber'})),double(dm1(:,{xI,yI,'ImageNumber'})),'rows');
end

function [dm,dropI]=getSlideOutlinesPerImage(dm, xI, yI)
RES=400;
xy=nan(2,size(dm,1));
xy(1,:)=round( RES*double(dm(:,xI))/nanmax(dm(:,xI)) )';
xy(2,:)=round( RES*double(dm(:,yI))/nanmax(dm(:,yI)) )';

keepI=[]; dropI=[];
for index=[1,2;2,1]
    for x = unique(xy(index(1),:))
        ii=find(x==xy(index(1),:));
        ulQ=quantile(xy(index(2),ii),[0.05,0.95]);
        ij=find(xy(index(2),:)<=ulQ(1) | xy(index(2),:)>=ulQ(2));
        keepI=[keepI,intersect(ii,ij)];
        %%now get outliers
        ulQ=quantile(xy(index(2),ii),[0.0001,0.9999]);
        ij=find(xy(index(2),:)<=ulQ(1) | xy(index(2),:)>=ulQ(2));
        dropI=[dropI,intersect(ii,ij)];
    end
end
keepI=setdiff(keepI,dropI);
dm=dm(unique(keepI),:);
end