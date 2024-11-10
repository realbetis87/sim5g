TModel=app.TModel;
boundary=TModel.Boundaries(app.BoundaryIndices);
Nv=numel(boundary.Vertices);Ne=numel(boundary.Facets);
p=zeros(2,Nv);t=zeros(4,Ne);
for ii=1:Nv,vertex=TModel.Vertices(boundary.Vertices(ii));
    p(1,vertex.Index2D)=vertex.Y;p(2,vertex.Index2D)=vertex.Z;
end
for ii=1:Ne,facet=TModel.Facets(boundary.Facets(ii));
    t(1,ii)=facet.Vertices2D(1);t(2,ii)=facet.Vertices2D(2);t(3,ii)=facet.Vertices2D(3);t(4,ii)=1;
end
save("2DMesh","t","p");
