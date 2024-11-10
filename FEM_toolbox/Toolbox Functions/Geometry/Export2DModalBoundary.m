function [TModel] = Export2DModalBoundary(TModel,BoundaryIdentifier),boundary=TModel.Boundaries(BoundaryIdentifier);
    for ii=1:numel(boundary.Vertices),TModel.Vertices(boundary.Vertices(ii)).Index2D=ii;end
    NodeToNode=zeros(numel(boundary.Vertices),numel(boundary.Vertices));
    for ii=1:numel(boundary.Edges),edge=TModel.Edges(boundary.Edges(ii));edge.Index2D=ii;edge.Vertices2D=[TModel.Vertices(edge.Vertices).Index2D];TModel.Edges(boundary.Edges(ii))=edge;
        NodeToNode(edge.Vertices2D(1),edge.Vertices2D(2))=ii;NodeToNode(edge.Vertices2D(2),edge.Vertices2D(1))=-ii;TModel.Edges(boundary.Edges(ii))=edge;
    end
    for ii=1:numel(boundary.Facets),facet=TModel.Facets(boundary.Facets(ii));facet.Index2D=ii;facet.Vertices2D=[TModel.Vertices(facet.Vertices).Index2D];
        facet.Medium2D=Medium();facet.Medium2D=TModel.Domains(TModel.Elements(facet.InElement(1)).SubDomain).Medium;TModel.Facets(boundary.Facets(ii))=facet;TModel.Facets(ii)=facet;
    end
end


