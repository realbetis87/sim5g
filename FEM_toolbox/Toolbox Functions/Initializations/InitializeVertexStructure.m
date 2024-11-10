function [Structure] = InitializeVertexStructure(Structure),BoundaryVertexCounter=0;
    Structure.NumberOfVertices=size(Structure.model.Mesh.Nodes,2);Structure.Vertices=Vertex.empty(0,Structure.NumberOfVertices);
    LV=zeros(Structure.NumberOfLineBoundaries,Structure.NumberOfVertices);BV=zeros(Structure.NumberOfBoundaries,Structure.NumberOfVertices);DV=zeros(Structure.NumberOfDomains,Structure.NumberOfVertices);
    for ii=1:Structure.NumberOfLineBoundaries,nodes=findNodes(Structure.model.Mesh,"region","Edge",ii);LV(ii,nodes)=1;Structure.LineBoundaries(ii).Vertices=nodes;end
    for ii=1:Structure.NumberOfBoundaries,nodes=findNodes(Structure.model.Mesh,"region","Face",ii);BV(ii,nodes)=1;Structure.Boundaries(ii).Vertices=nodes;end
    for ii=1:Structure.NumberOfDomains,nodes=findNodes(Structure.model.Mesh,"region","Cell",ii);DV(ii,nodes)=1;Structure.Domains(ii).Vertices=nodes;end
    for ii=1:Structure.NumberOfVertices,inDomain=find(DV(:,ii));onLine=find(LV(:,ii));onBoundary=find(BV(:,ii));
        if(isempty(onLine) && ~isempty(onBoundary)),Structure.Vertices(ii)=Vertex(Structure,ii,Structure.model.Mesh.Nodes(1,ii),Structure.model.Mesh.Nodes(2,ii),Structure.model.Mesh.Nodes(3,ii),inDomain,onBoundary);
        elseif(isempty(onLine) && isempty(onBoundary)),Structure.Vertices(ii)=Vertex(ii,Structure.model.Mesh.Nodes(1,ii),Structure.model.Mesh.Nodes(2,ii),Structure.model.Mesh.Nodes(3,ii),inDomain);
        else,Structure.Vertices(ii)=Vertex(Structure,ii,Structure.model.Mesh.Nodes(1,ii),Structure.model.Mesh.Nodes(2,ii),Structure.model.Mesh.Nodes(3,ii),inDomain,onBoundary,onLine);
        end
        if(Structure.Vertices(ii).OnExterior),BoundaryVertexCounter=BoundaryVertexCounter+1;end
    end
    for ii=1:Structure.NumberOfElements,element=Structure.Elements(ii);Structure.Elements(ii)=Structure.Elements(ii).UpdateElement(Structure.Vertices);
        for kk=1:4,vertex=Structure.Vertices(element.Vertices(kk));vertex.InElement(end+1)=element.Index;Structure.Vertices(element.Vertices(kk))=vertex;end
    end,Structure.NumberOfBoundaryVertices=BoundaryVertexCounter;
    [Structure] = DetermineBoundaryGeometry(Structure);
end