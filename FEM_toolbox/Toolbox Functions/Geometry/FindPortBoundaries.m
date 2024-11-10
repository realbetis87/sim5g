function [TModel,BoundaryIndices] = FindPortBoundaries(TModel,BoundaryIndex),boundary=TModel.Boundaries(BoundaryIndex);error=1000*eps;
    if(boundary.Type=="POR")
            d1=[TModel.Boundaries.Type]=="POR";d2=[TModel.Boundaries.Axis]==boundary.Axis;d3=zeros(numel(TModel.Boundaries),1);
            for ii=1:numel(TModel.Boundaries),if(~isempty(TModel.Boundaries(ii).Position)),if(abs(TModel.Boundaries(ii).Position-boundary.Position)<=error),d3(ii)=1;end,end
            end,d=d1 & d2 & d3';boundaryIndices=find(d);BoundaryIndices=boundaryIndices;
    elseif(boundary.Type=="DIR")
             d1=[TModel.Boundaries.Type]=="DIR";d2=[TModel.Boundaries.Axis]==boundary.Axis;d3=zeros(numel(TModel.Boundaries),1);
            for ii=1:numel(TModel.Boundaries),if(~isempty(TModel.Boundaries(ii).Position)),if(abs(TModel.Boundaries(ii).Position-boundary.Position)<=error),d3(ii)=1;end,end
            end,d=d1 & d2 & d3';boundaryIndices=find(d);BoundaryIndices=boundaryIndices;
    end,TModel=Prepare2DBoundaries(TModel,boundaryIndices);
end
function TModel = Prepare2DBoundaries(TModel,boundaryIndices),vc=0;      
    for ii=1:numel(boundaryIndices),boundary=TModel.Boundaries(boundaryIndices(ii));
        for jj=1:numel(boundary.Vertices),vertex=TModel.Vertices(boundary.Vertices(jj));
            if(isempty(vertex.Index2D)),vc=vc+1;vertex.Index2D=vc;end,TModel.Vertices(boundary.Vertices(jj))=vertex;
        end
        for jj=1:numel(boundary.Facets),facet=TModel.Facets(boundary.Facets(jj));
            if(isscalar(facet.InElement)),element=TModel.Elements(facet.InElement);else,element=TModel.Elements(facet.InElement(1));end,facet.Medium2D=TModel.Domains(element.SubDomain).Medium;
            verts=[TModel.Vertices(facet.Vertices)];vf1=verts(1).Index;vf2=verts(2).Index;vf3=verts(3).Index;
            facet.Vertices2D(1)=verts(1).Index2D;facet.Vertices2D(2)=verts(2).Index2D;facet.Vertices2D(3)=verts(3).Index2D;
            for kk=1:6,edge=TModel.Edges(element.Edges(kk));ve1=edge.Vertices(1);ve2=edge.Vertices(2);
                if(ve1==vf1 && ve2==vf2),facet.Edges(1)=edge.Index;facet.EdgeSigns(1)=1;
                elseif(ve1==vf2 && ve2==vf1),facet.Edges(1)=edge.Index;facet.EdgeSigns(1)=-1;
                elseif(ve1==vf2 && ve2==vf3),facet.Edges(2)=edge.Index;facet.EdgeSigns(2)=1;
                elseif(ve1==vf3 && ve2==vf2),facet.Edges(2)=edge.Index;facet.EdgeSigns(2)=-1;
                elseif(ve1==vf3 && ve2==vf1),facet.Edges(3)=edge.Index;facet.EdgeSigns(3)=1;
                elseif(ve1==vf1 && ve2==vf3),facet.Edges(3)=edge.Index;facet.EdgeSigns(3)=-1;
                end
            end,TModel.Facets(boundary.Facets(jj))=facet;
        end
    end
end