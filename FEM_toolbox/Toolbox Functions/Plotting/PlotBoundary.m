function [] = PlotBoundary(Structure,Type,Identifier)
    if(Type=="Vertices"),PlotBoundaryVertices(Structure,Identifier);
    elseif(Type=="Edges"),PlotBoundaryEdges(Structure,Identifier);
    elseif(Type=="Facets"),PlotBoundaryFacets(Structure,Identifier);
    end
end

function []  = PlotBoundaryVertices(Structure,Identifier)
    boundary=Structure.Boundaries(Identifier);
    for ii=1:numel(boundary.Vertices),vertex=Structure.Vertices(boundary.Vertices(ii));plot3(vertex.X,vertex.Y,vertex.Z,'o','Color','k');hold on;end
end
function [] = PlotBoundaryEdges(Structure,Identifier)
    boundary=Structure.Boundaries(Identifier);
    for ii=1:numel(boundary.Edges),edge=Structure.Edges(boundary.Edges(ii));vertex1=Structure.Vertices(edge.Vertices(1));vertex2=Structure.Vertices(edge.Vertices(2));plot3([vertex1.X vertex2.X],[vertex1.Y vertex2.Y],[vertex1.Z vertex2.Z],'Color','r');hold on;end
end
function [] = PlotBoundaryFacets(Structure,Identifier)
    boundary=Structure.Boundaries(Identifier);
    for ii=1:numel(boundary.Facets),facet=Structure.Facets(boundary.Facets(ii));vertex1=Structure.Vertices(facet.Vertices(1));vertex2=Structure.Vertices(facet.Vertices(2));vertex3=Structure.Vertices(facet.Vertices(3));
        P1=[vertex1.X,vertex2.X,vertex3.X];P2=[vertex1.Y,vertex2.Y,vertex3.Y];P3=[vertex1.Z,vertex2.Z,vertex3.Z];fill3(P1,P2,P3,[0.8500 0.3250 0.0980],'FaceAlpha',0.4);hold on;
        plot3([vertex1.X vertex2.X],[vertex1.Y vertex2.Y],[vertex1.Z vertex2.Z],'Color','k');hold on;
        plot3([vertex2.X vertex3.X],[vertex2.Y vertex3.Y],[vertex2.Z vertex3.Z],'Color','k');hold on;
        plot3([vertex3.X vertex1.X],[vertex3.Y vertex1.Y],[vertex3.Z vertex1.Z],'Color','k');hold on;
    end
end