function dmnew=addColumn(dmobj,colname, initval)
import bioma.data.*

dmnew=dmobj;
colname=char(colname);
colsnew =colnames(dmobj);
colsnew{length(colsnew)+1}=colname;
vals=[double(dmobj),ones(size(dmobj,1),1)*initval];
if ~any(strcmp(colnames(dmobj),colname))
    dmnew=DataMatrix(vals,rownames(dmobj),colsnew);
end

end