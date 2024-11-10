function [toolboxModel] = InitializeFacetStructure(toolboxModel,NodeToNode)
    toolboxModel.NumberOfFacets=toolboxModel.NumberOfEdges-toolboxModel.NumberOfVertices+toolboxModel.NumberOfElements-1;
    toolboxModel.Facets=Facet.empty(0,toolboxModel.NumberOfFacets);FirstNode=[1 2 3 4];SecondNode=[2 4 4 2];ThirdNode=[3 3 1 1];FacetCounter=0;
    for ie=1:toolboxModel.NumberOfElements
        for ii=1:4,Node1=toolboxModel.Elements(ie).Vertices(FirstNode(ii));Node2=toolboxModel.Elements(ie).Vertices(SecondNode(ii));Node3=toolboxModel.Elements(ie).Vertices(ThirdNode(ii));
             if(toolboxModel.Elements(ie).Facets(ii)==0),FacetCounter=FacetCounter+1;toolboxModel.Facets(FacetCounter)=Facet(toolboxModel,FacetCounter,Node1,Node2,Node3,ie);toolboxModel=toolboxModel.Facets(FacetCounter).UpdateStructureFacet(toolboxModel);toolboxModel.Elements(ie).Facets(ii)=FacetCounter;toolboxModel.Elements(ie).FacetSigns(ii)=1;
             end
             %if(FacetCounter==toolboxModel.NumberOfFacets),break;end
        end
    end
    toolboxModel.NumberOfFacets=FacetCounter;
    for ii=1:toolboxModel.NumberOfFacets,facet=toolboxModel.Facets(ii);
        facet.Edges(1)=abs(NodeToNode(facet.Vertices(1),facet.Vertices(2)));facet.EdgeSigns(1)=sign(NodeToNode(facet.Vertices(1),facet.Vertices(2)));
        facet.Edges(2)=abs(NodeToNode(facet.Vertices(2),facet.Vertices(3)));facet.EdgeSigns(2)=sign(NodeToNode(facet.Vertices(2),facet.Vertices(3)));
        facet.Edges(3)=abs(NodeToNode(facet.Vertices(3),facet.Vertices(1)));facet.EdgeSigns(3)=sign(NodeToNode(facet.Vertices(3),facet.Vertices(1)));
        toolboxModel.Facets(ii)=facet;toolboxModel=toolboxModel.Facets(ii).UpdateFacetBoundary(toolboxModel);
    end
    for ii=1:numel(toolboxModel.Boundaries),boundary=toolboxModel.Boundaries(ii);
        for jj=1:numel(boundary.Facets),toolboxModel.Facets(boundary.Facets(jj)).OnBoundary=ii;end
        for jj=1:numel(boundary.Edges),if(numel(toolboxModel.Edges(boundary.Edges(jj)).OnBoundary)>1),toolboxModel.Edges(boundary.Edges(jj)).OnBoundary(end+1)=ii;else,toolboxModel.Edges(boundary.Edges(jj)).OnBoundary=ii;end,end
    end
    for ii=1:numel(toolboxModel.Facets)
    end
end
   


