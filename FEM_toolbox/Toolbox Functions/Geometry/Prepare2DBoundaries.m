function [TModel,BoundaryIndices,BoundaryExcitation,LineBoundaryIndices] = Prepare2DBoundaries(TModel,BoundaryIndex),vc=0;LineBoundaryIndices=[];   
         [TModel,BoundaryIndices]=FindPortBoundaries(TModel,BoundaryIndex);
         vertices=[TModel.Boundaries(BoundaryIndices(1)).Vertices];
         edges=[TModel.Boundaries(BoundaryIndices(1)).Edges];
         facets=[TModel.Boundaries(BoundaryIndices(1)).Facets];
         for ii=1:numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(ii));
             temp=vertices;vertices=zeros(numel(temp)+numel(boundary.Vertices),1);vertices(1:numel(temp))=temp;vertices(numel(temp)+1:end)=boundary.Vertices;
             temp=edges;edges=zeros(numel(temp)+numel(boundary.Edges),1);edges(1:numel(temp))=temp;edges(numel(temp)+1:end)=boundary.Edges;
             temp=facets;facets=zeros(numel(temp)+numel(boundary.Facets),1);facets(1:numel(temp))=temp;facets(numel(temp)+1:end)=boundary.Facets;
             for jj=1:numel(boundary.Lines)
                    if(nnz(intersect(LineBoundaryIndices,boundary.Lines(jj)))==0),temp=LineBoundaryIndices;LineBoundaryIndices=zeros(numel(temp)+1,1);LineBoundaryIndices(1:end-1)=temp;LineBoundaryIndices(end)=boundary.Lines(jj);end
             end
         end
         vertices=unique(vertices);edges=unique(edges);facets=unique(facets);
         for ii=1:numel(vertices),TModel.Vertices(vertices(ii)).Index2D=ii;end
         for ii=1:numel(edges),TModel.Edges(edges(ii)).Index2D=ii;end
         for ii=1:numel(facets),facet=TModel.Facets(facets(ii));
             if(isscalar(facet.InElement)),element=TModel.Elements(facet.InElement);
             else,element=TModel.Elements(facet.InElement(1));
             end,facet.Medium2D=TModel.Domains(element.SubDomain).Medium;
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
            end,TModel.Facets(facets(ii))=facet;
         end
    Type=boundary.Type;BoundaryExcitation=Excitation(BoundaryIndices,Type);
    BoundaryExcitation.LineBoundaries=LineBoundaryIndices;
    BoundaryExcitation.BoundaryIndices=BoundaryIndices;
    BoundaryExcitation.Vertices=vertices;
    BoundaryExcitation.Edges=edges;
    BoundaryExcitation.Facets=facets;
    BoundaryExcitation=Determine2DModal_Equation(TModel,BoundaryExcitation,BoundaryIndices);
end
%{
function [TModel,BoundaryIndices,BoundaryExcitation,LineBoundaryIndices] = Prepare2DBoundaries(TModel,BoundaryIndex),vc=0;LineBoundaryIndices=[];   
         [TModel,BoundaryIndices]=FindPortBoundaries(TModel,BoundaryIndex);
%---------------- 2D Vertex Index Enumeration -----------------------------
    for ii=1:numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(ii));
        for jj=1:numel(boundary.Vertices),vertex=TModel.Vertices(boundary.Vertices(jj));
            if(isempty(vertex.Index2D)),vc=vc+1;vertex.Index2D=vc;end,TModel.Vertices(boundary.Vertices(jj))=vertex;
        end
        %------------------- 2D Edge Matching -----------------------------
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
        %------------- Line Boundary Identification -----------------------
        for jj=1:numel(boundary.Lines)
            if(nnz(intersect(LineBoundaryIndices,boundary.Lines(jj)))==0),temp=LineBoundaryIndices;LineBoundaryIndices=zeros(numel(temp)+1,1);LineBoundaryIndices(1:end-1)=temp;LineBoundaryIndices(end)=boundary.Lines(jj);end
        end
    end,Type=boundary.Type;BoundaryExcitation=Excitation(BoundaryIndices,Type);BoundaryExcitation.LineBoundaries=LineBoundaryIndices;BoundaryExcitation.BoundaryIndices=BoundaryIndices;
    Edges=zeros(numel(TModel.Edges),1);Facets=zeros(numel(TModel.Facets),1);Vertices=zeros(numel(TModel.Vertices),1);
    for ib=1:numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(ib));
        for ii=1:numel(boundary.Vertices),vertex=TModel.Vertices(boundary.Vertices(ii));if(Vertices(vertex.Index)==0),Vertices(vertex.Index)=1;BoundaryExcitation=BoundaryExcitation.AddVertex(vertex.Index);TModel.Vertices(vertex.Index).Index2D=ii;end,end
        for ii=1:numel(boundary.Edges),edge=TModel.Edges(boundary.Edges(ii));if(Edges(edge.Index)==0),Edges(edge.Index)=1;BoundaryExcitation=BoundaryExcitation.AddEdge(edge.Index);TModel.Edges(edge.Index).Index2D=ii;end,end
        for ii=1:numel(boundary.Facets),facet=TModel.Facets(boundary.Facets(ii));if(Facets(facet.Index)==0),Facets(facet.Index)=1;BoundaryExcitation=BoundaryExcitation.AddFacet(facet.Index);TModel.Facets(facet.Index).Index2D=ii;end,end
    end
    BoundaryExcitation=Determine2DModal_Equation(TModel,BoundaryExcitation,BoundaryIndices);
end
%}
function [BoundaryExcitation] = Determine2DModal_Equation(TModel,BoundaryExcitation,BoundaryIndices),bianFlag=false;orthFlag=false;
        for ii=1:numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(ii));element=TModel.Facets(boundary.Facets(1));medium=element.Medium2D;
                if(medium.Type~="Iso"),if(medium.Type=="Bian"),bianFlag=true;
                                       elseif(CheckNonOrthotropicMedium(boundary,medium)),orthFlag=true;
                                       end
                end
        end
        if(bianFlag==true),Equation=3;
        elseif(bianFlag==false && orthFlag==true),Equation=2;
        else,Equation=1;
        end
        BoundaryExcitation.Assembly=FEMAssembly("EigenMode","2D",Equation);
end
%------------------- Orthotropic Check ------------------------------------
function [res] = CheckNonOrthotropicMedium(boundary,medium),res=false;
    if(medium.IsDispersive)
        for ii=1:medium.FRange.NF,r1=CheckNonOrthotropicTensor(boundary,medium.Epsilon{ii});r2=CheckNonOrthotropicTensor(boundary,medium.Mu{ii});
            if(r1 || r2),res=true;end
        end
    else,r1=CheckNonOrthotropicTensor(boundary,medium.Epsilon);r2=CheckNonOrthotropicTensor(boundary,medium.Mu);if(r1 || r2),res=true;end
    end
end
function [res] = CheckNonOrthotropicTensor(boundary,tensor),res=false;
    switch boundary.Axis
        case  1,if(tensor(1,2)~=0 || tensor(1,3)~=0 || tensor(3,1)~=0 || tensor(2,1)~=0),res=true;end
        case -1,if(tensor(1,2)~=0 || tensor(1,3)~=0 || tensor(3,1)~=0 || tensor(2,1)~=0),res=true;end
        case  2,if(tensor(2,1)~=0 || tensor(2,3)~=0 || tensor(1,2)~=0 || tensor(3,2)~=0),res=true;end
        case -2,if(tensor(2,1)~=0 || tensor(2,3)~=0 || tensor(1,2)~=0 || tensor(3,2)~=0),res=true;end
        case  3,if(tensor(3,1)~=0 || tensor(3,2)~=0 || tensor(2,3)~=0 || tensor(1,3)~=0),res=true;end
        case -3,if(tensor(3,1)~=0 || tensor(3,2)~=0 || tensor(2,3)~=0 || tensor(1,3)~=0),res=true;end
    end
end