function [] = UpdateFacets_GUI(app,facet)
    Vertex1=app.Vertices(facet.Vertices(1));Vertex2=app.Vertices(facet.Vertices(2));Vertex3=app.Vertices(facet.Vertices(3));
    CommonElements=intersect(Vertex1.InElement,Vertex2.InElement);CommonElements=intersect(CommonElements,Vertex3.InElement);facet.InElement=CommonElements;
    for ii=1:numel(CommonElements),element=app.Elements(CommonElements(ii));[facetSign,facetPosition]=FacetPositions(element,facet);element.Facets(facetPosition)=facet.Index;element.FacetSigns(facetPosition)=facetSign;app.Elements(CommonElements(ii))=element;end
    facet=facet.UpdateFacet(app.Vertices);facet=CheckBoundary(app,facet);app.Facets(facet.Index)=facet;
           

end
function [Sign,Position] = FacetPositions(element,facet),Sign=1;Position=1;V=facet.Vertices;FacetPositions=[1 2 3 4;2 4 4 2;3 3 1 1;];ElementVertices=element.Vertices;
    for ii=1:4,ElementV=ElementVertices(FacetPositions(:,ii));Check=double(eq(ElementV,V));if(nnz(Check)==3),Position=ii;Permutations=nnz(eye(3)-Check);if(nnz(Permutations)==0 || nnz(Permutations)==6),Sign=1;else,Sign=-1;end,end,end
end

function [facet] = CheckBoundary(app,facet),V1=app.Vertices(facet.Vertices(1));V2=app.Vertices(facet.Vertices(2));V3=app.Vertices(facet.Vertices(3));
       if(~isempty(V1.OnBoundary) && ~isempty(V2.OnBoundary) && ~isempty(V1.OnBoundary)),res1=intersect(V1.OnBoundary,V2.OnBoundary);
          if(~isempty(res1) && ~isempty(intersect(res1,V3.OnBoundary))),facet.OnBoundary=min(intersect(res1,V3.OnBoundary));
              %app.Boundaries(facet.OnBoundary)=app.Boundaries(facet.OnBoundary).AddFacet(facet.Index);
          end
       end
       
end
