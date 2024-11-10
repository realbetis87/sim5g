rr=1:546; boundary=TModel.Boundaries(5);

for ii=1:numel(boundary.Edges)
edge=TModel.Edges(boundary.Edges(ii));
v1=TModel.Vertices(edge.Vertices(1));v2=TModel.Vertices(edge.Vertices(2));
for jj=1:221
if((v1.Index2D == Edge_Vertices(jj,1) && v2.Index2D == Edge_Vertices(jj,2)) || (v2.Index2D == Edge_Vertices(jj,1) && v1.Index2D == Edge_Vertices(jj,2)))
if(edge.IndexH~=0 && Index_H(jj)~=0),rr(Index_H(jj))=edge.IndexH;
    disp(1);
end
if(edge.IndexE~=0 && Index_E(jj)~=0),rr(Index_E(jj))=edge.IndexE;
    disp(2);
end
end
end
end