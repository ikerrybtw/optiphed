function dmO=removeNoiseFromImage(dm, xI, yI)
dmO=[];
if any(strcmp(colnames(dm),'ImageNumber'))
    for img=unique(double(dm(:,'ImageNumber')))'
        dm1=removeNoisePerImage(dm(dm(:,'ImageNumber')==img,:),xI,yI);
        dmO=[dmO;dm1];
    end
else
    dmO=removeNoisePerImage(dm,xI,yI);
end
end

function dm=removeNoisePerImage(dm, xI, yI)
RES=400;
xy=nan(2,size(dm,1));
xy(1,:)=round( RES*double(dm(:,xI))/nanmax(dm(:,xI)) )';
xy(2,:)=round( RES*double(dm(:,yI))/nanmax(dm(:,yI)) )';

dropI=[];
for index=[1,2;2,1]
    uniqueAx=unique(xy(index(1),:));
    counts=grpstats(xy(index(1),:),xy(index(1),:),'numel');
    for x = uniqueAx
        ii=find(x==xy(index(1),:));
        if length(ii)<0.01*max(counts); %0.1*mean(counts)
            dropI=[dropI,ii];
        end
    end
end
dropI=unique(dropI);
if ~isempty(dropI)
    disp(['Removing noisy cell-artifacts: ',num2str(length(dropI)),' dots (',num2str(100*length(dropI)/size(dm,1)),'%)'])
    dm(dropI,:)=[];
end
end