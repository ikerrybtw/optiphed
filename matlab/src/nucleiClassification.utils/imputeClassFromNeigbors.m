function cl=imputeClassFromNeigbors(hande_dm,maxD,xI,yI,cl)
nans=sum(isnan(cl));
last=length(cl);
while nans<last
    for i=find(isnan(cl))'
        d=pdist2(double(hande_dm(i,{xI,yI})),double(hande_dm(:,{xI,yI})));
        [a,ia]=sort(d);
        ia=ia(a<maxD);%%k closest neighbors
        cl(i)=round(nanmean(cl(ia)));
    end
    last=nans;
    nans=sum(isnan(cl));
end
end