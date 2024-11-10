function [Structure,NodeToNode] = InitializeEdgeStructure(Structure)
    Structure.NumberOfEdges=Structure.NumberOfElements+Structure.NumberOfVertices+Structure.NumberOfBoundaryVertices-3;Structure.Edges=Edge.empty(0,Structure.NumberOfEdges);
    NodeToNode=spalloc(Structure.NumberOfVertices,Structure.NumberOfVertices,6*Structure.NumberOfVertices);FirstNode=[1 1 1 2 2 3];SecondNode=[2 3 4 3 4 4];EdgeCounter=0;
    for ie=1:Structure.NumberOfElements,element=Structure.Elements(ie);
        for in=1:6,i1=FirstNode(in);i2=SecondNode(in);Node1=element.Vertices(i1);Node2=element.Vertices(i2);
            if(NodeToNode(Node1,Node2)==0),EdgeCounter=EdgeCounter+1;NodeToNode(Node1,Node2)=EdgeCounter;NodeToNode(Node2,Node1)=-EdgeCounter;
                Structure.Edges(EdgeCounter)=Edge(Structure,EdgeCounter,Node1,Node2,ie);
                Structure.Elements(ie)=Structure.Elements(ie).AddEdge(in,EdgeCounter,1);
                Structure=BoundaryCheck(Structure.Edges(EdgeCounter),Structure);
            else
                Structure.Edges(abs(NodeToNode(Node1,Node2))).InElement(end+1)=ie;
                Structure.Elements(ie)=Structure.Elements(ie).AddEdge(in,abs(NodeToNode(Node1,Node2)),sign(NodeToNode(Node1,Node2)));
            end
        end
    end
end

function [Structure] = BoundaryCheck(edge,Structure)
    if(~isempty(edge.OnLine) && edge.OnLine~=0),Structure.LineBoundaries(edge.OnLine)=Structure.LineBoundaries(edge.OnLine).AddEdge(edge.Index);end
    if(~isempty(edge.OnBoundary) && all(edge.OnBoundary~=0)),for ee=1:numel(edge.OnBoundary),bi=edge.OnBoundary(ee);Structure.Boundaries(bi)=Structure.Boundaries(bi).AddEdge(edge.Index);end
    end
end