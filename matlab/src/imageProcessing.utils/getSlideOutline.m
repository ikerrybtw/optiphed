function keepI=getSlideOutline(dm, xI, yI,side,thickness)
if ~exist('thickness','var')
    thickness=0.05;
end
if  strcmp(side,'all')
    keepI=getSlideOutlinesPerImage(dm, xI, yI, thickness);
elseif strcmp(side,'upper')
    keepI=getSlideOutlineOnSide(dm, xI, yI,[1,2],thickness);
elseif strcmp(side,'lower')
    keepI=getSlideOutlineOnSide(-dm, xI, yI,[1,2],thickness);
elseif strcmp(side,'right')
    keepI=getSlideOutlineOnSide(dm, xI, yI,[2,1],thickness);
elseif strcmp(side,'left')
    keepI=getSlideOutlineOnSide(-dm, xI, yI,[2,1],thickness);
end
end


function keepI=getSlideOutlineOnSide(dm, xI, yI,index,thickness)
RES=400;
xy=nan(2,size(dm,1));
xy(1,:)=round( RES*double(dm(:,xI))/nanmax(dm(:,xI)) )';
xy(2,:)=round( RES*double(dm(:,yI))/nanmax(dm(:,yI)) )';

keepI=[];
% index=[1,2];
for x = unique(xy(index(1),:))
    ii=find(x==xy(index(1),:));
    ulQ=quantile(xy(index(2),ii),1-thickness);
    ij=find(xy(index(2),:)>=ulQ);
    keepI=[keepI,intersect(ii,ij)];
end
keepI=unique(keepI,'stable');
end

function keepI=getSlideOutlinesPerImage(dm, xI, yI, thickness)
RES=400;
xy=nan(2,size(dm,1));
xy(1,:)=round( RES*double(dm(:,xI))/nanmax(dm(:,xI)) )';
xy(2,:)=round( RES*double(dm(:,yI))/nanmax(dm(:,yI)) )';

keepI=[];
for index=[1,2;2,1]
    for x = unique(xy(index(1),:))
        ii=find(x==xy(index(1),:));
        ulQ=quantile(xy(index(2),ii),[thickness,1-thickness]);
        ij=find(xy(index(2),:)<=ulQ(1) | xy(index(2),:)>=ulQ(2));
        keepI=[keepI,intersect(ii,ij)];
    end
end
keepI=unique(keepI,'stable');
end