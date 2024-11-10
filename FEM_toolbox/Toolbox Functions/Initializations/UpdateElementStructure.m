function [toolboxModel] = UpdateElementStructure(toolboxModel)
    ve1=[1 1 1 2 2 3];ve2=[2 3 4 3 4 4];
    for ie=1:numel(toolboxModel.Elements),element=toolboxModel.Elements(ie);element.FacetSigns=zeros(4,1);
        for ii=1:6,edge=toolboxModel.Edges(element.Edges(ii));
           if(element.Vertices(ve1(ii))==edge.Vertices(1) && element.Vertices(ve2(ii))==edge.Vertices(2)),element.EdgeSigns(ii)=1;
           elseif(element.Vertices(ve2(ii))==edge.Vertices(1) && element.Vertices(ve1(ii))==edge.Vertices(2)),element.EdgeSigns(ii)=-1;
           end
        end
        for jj=1:4,facet=toolboxModel.Facets(element.Facets(jj));v1=toolboxModel.Vertices(facet.Vertices(1));v2=toolboxModel.Vertices(facet.Vertices(2));v3=toolboxModel.Vertices(facet.Vertices(3));
            
        end
    end
end

function sign = DetermineFacetSing(element,facet,facetIndex),v1=element.Vertices(1);v2=element.Vertices(2);v3=element.Vertices(3);v4=element.Vertices(4);
    switch facetIndex
        case 1
        case 2
        case 3
        case 4
    end
end

function [Sign,Position] = FacetPositionSign(element,Vertex1,Vertex2,Vertex3),Sign=0;Position=1;V=[Vertex1.Index Vertex2.Index Vertex3.Index];FacetPositions=[1 2 3 4;2 4 4 2;3 3 1 1;];
     ElementVertices=element.Vertices;
    for ii=1:4,ElementV=ElementVertices(FacetPositions(:,ii));Check=double(eq(ElementV,V));if(nnz(Check)==3),Position=ii;Permutations=nnz(eye(3)-Check);if(nnz(Permutations)==0 || nnz(Permutations)==6),Sign=1;else,Sign=-1;end,end,end
end