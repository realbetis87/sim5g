function [TModel] = EnumerateDoFs_EB(TModel),UknownCounter=0;KnownCounter=0;
%===================== Edge DoF Enumeration ===============================
%-------------------- Periodic Boundaries ---------------------------------
    for ib=1:numel(TModel.Boundaries),boundary=TModel.Boundaries(ib);
        if(boundary.Type=="PBC")
            if(boundary.Master==1)
                for ie=1:numel(boundary.Edges),edge=TModel.Edges(boundary.Edges(ie));
                      if(isempty(edge.UknownIndex)),UknownCounter=UknownCounter+1;
                            edge.UknownIndex=UknownCounter;edge.KnownIndex=0;
                            slaveIndex=findEdgePeriodicPair(TModel,edge,boundary);edge.PPair=abs(slaveIndex);slave_edge=TModel.Edges(edge.PPair);slave_edge.KnownIndex=0;
                            slave_edge.UknownIndex=sign(slaveIndex)*UknownCounter;slave_edge.PPair=edge.Index;TModel.Edges(edge.PPair)=slave_edge;TModel.Edges(boundary.Edges(ie))=edge;
                       end
                end
             end
        end
    end
%-------------------- Boundaries with DoFs --------------------------------
    for ib=1:numel(TModel.Boundaries),boundary=TModel.Boundaries(ib);
        if(boundary.Id==0)
            for ie=1:numel(boundary.Edges),edge=TModel.Edges(boundary.Edges(ie));
                if(isempty(edge.UknownIndex)),UknownCounter=UknownCounter+1;edge.UknownIndex=UknownCounter;edge.KnownIndex=0;TModel.Edges(boundary.Edges(ie))=edge;end
            end
        elseif(boundary.Id~=0 && boundary.Type=="DIR")
            for ie=1:numel(boundary.Edges),edge=TModel.Edges(boundary.Edges(ie));
                if(isempty(edge.UknownIndex)),KnownCounter=KnownCounter+1;edge.KnownIndex=KnownCounter;edge.UknownIndex=0;TModel.Edges(boundary.Edges(ie))=edge;end
            end
        end
    end
%------------------ Boundaries without DoFs -------------------------------
    for ib=1:numel(TModel.Boundaries),boundary=TModel.Boundaries(ib);
        if(boundary.Id~=0 && boundary.Type~="DIR")
            for ie=1:numel(boundary.Edges),edge=TModel.Edges(boundary.Edges(ie));
                if(isempty(edge.UknownIndex)),edge.KnownIndex=0;edge.UknownIndex=0;TModel.Edges(boundary.Edges(ie))=edge;end
            end
        end
    end
%-----------------------------  Rest  -------------------------------------
    for ie=1:numel(TModel.Edges),edge=TModel.Edges(ie);
        if(isempty(edge.UknownIndex)),elements=[TModel.Elements(edge.InElement)];subDoms=[elements.SubDomain];
            if(any([TModel.Domains(subDoms).IsEmpty]==true)),edge.KnownIndex=0;edge.UknownIndex=0;TModel.Edges(ie)=edge;
            else,edge.KnownIndex=0;UknownCounter=UknownCounter+1;edge.UknownIndex=UknownCounter;TModel.Edges(ie)=edge;
            end
        end
    end
%===================== Facet DoF Enumeration ==============================
%-------------------- Periodic Boundaries ---------------------------------
    for ib=1:numel(TModel.Boundaries),boundary=TModel.Boundaries(ib);
        if(boundary.Type=="PBC")
            if(boundary.Master==1)
                for ie=1:numel(boundary.Facets),facet=TModel.Facets(boundary.Facets(ie));
                      if(isempty(facet.UknownIndex)),UknownCounter=UknownCounter+1;facet.UknownIndex=UknownCounter;facet.KnownIndex=0;
                            slaveIndex=findFacetPeriodicPair(TModel,facet,boundary);facet.PPair=abs(slaveIndex);slave_facet=TModel.Facets(facet.PPair);slave_facet.KnownIndex=0;
                            slave_facet.UknownIndex=sign(slaveIndex)*UknownCounter;slave_facet.PPair=facet.Index;TModel.Facets(facet.PPair)=slave_facet;TModel.Facets(boundary.Facets(ie))=facet;
                       end
                end
             end
        end
    end
%-------------------- Boundaries with DoFs --------------------------------
    for ib=1:numel(TModel.Boundaries),boundary=TModel.Boundaries(ib);
        if(boundary.Id==0)
            for ie=1:numel(boundary.Facets),facet=TModel.Facets(boundary.Facets(ie));
                if(isempty(facet.UknownIndex)),UknownCounter=UknownCounter+1;facet.UknownIndex=UknownCounter;facet.KnownIndex=0;TModel.Facets(boundary.Facets(ie))=facet;end
            end
        elseif(boundary.Id~=0 && boundary.Type=="DIR")
            for ie=1:numel(boundary.Facets),facet=TModel.Facets(boundary.Facets(ie));
                if(isempty(facet.UknownIndex)),KnownCounter=KnownCounter+1;facet.KnownIndex=KnownCounter;facet.UknownIndex=0;TModel.Facets(boundary.Facets(ie))=facet;end
            end
        end
    end
 %------------------ Boundaries without DoFs -------------------------------
    for ib=1:numel(TModel.Boundaries),boundary=TModel.Boundaries(ib);
        if(boundary.Id~=0 && boundary.Type~="DIR")
            for ie=1:numel(boundary.Facets),facet=TModel.Facets(boundary.Facets(ie));
                if(isempty(facet.UknownIndex)),facet.KnownIndex=0;facet.UknownIndex=0;TModel.Facets(boundary.Facets(ie))=facet;end
            end
        end
    end
%-----------------------------  Rest  -------------------------------------
    for ie=1:numel(TModel.Facets),facet=TModel.Facets(ie);
        if(isempty(facet.UknownIndex)),elements=[TModel.Elements(facet.InElement)];subDoms=[elements.SubDomain];
            if(any([TModel.Domains(subDoms).IsEmpty]==true)),facet.KnownIndex=0;facet.UknownIndex=0;TModel.Facets(ie)=facet;
            else,UknownCounter=UknownCounter+1;facet.KnownIndex=0;facet.UknownIndex=UknownCounter;TModel.Facets(ie)=facet;
            end 
         end
    end
end
%--------------------------------------------------------------------------
function [pair_index] = findEdgePeriodicPair(TModel,edge,boundary),slave_Boundary=TModel.Boundaries(boundary.Param);error=1000*eps;
    v1=TModel.Vertices(edge.Vertices(1));v2=TModel.Vertices(edge.Vertices(2));
    switch abs(boundary.Axis)
        case 1,for ii=1:numel(slave_Boundary.Edges),slave=TModel.Edges(slave_Boundary.Edges(ii));
                    sv1=TModel.Vertices(slave.Vertices(1));sv2=TModel.Vertices(slave.Vertices(2));
                    if(abs(v1.Y-sv1.Y)<=error && abs(v1.Z-sv1.Z)<=error && abs(v2.Y-sv2.Y)<=error && abs(v2.Z-sv2.Z)<=error),pair_index=slave.Index;
                    elseif(abs(v1.Y-sv2.Y)<=error && abs(v1.Z-sv2.Z)<=error && abs(v2.Y-sv1.Y)<=error && abs(v2.Z-sv1.Z)<=error),pair_index=-slave.Index;
                    end
                end
        case 2,for ii=1:numel(slave_Boundary.Edges),slave=TModel.Edges(slave_Boundary.Edges(ii));
                    sv1=TModel.Vertices(slave.Vertices(1));sv2=TModel.Vertices(slave.Vertices(2));
                    if(abs(v1.X-sv1.X)<=error && abs(v1.Z-sv1.Z)<=error && abs(v2.X-sv2.X)<=error && abs(v2.Z-sv2.Z)<=error),pair_index=slave.Index;
                    elseif(abs(v1.X-sv2.X)<=error && abs(v1.Z-sv2.Z)<=error && abs(v2.X-sv1.X)<=error && abs(v2.Z-sv1.Z)<=error),pair_index=-slave.Index;
                    end
                end
        case 3,for ii=1:numel(slave_Boundary.Edges),slave=TModel.Edges(slave_Boundary.Edges(ii));
                    sv1=TModel.Vertices(slave.Vertices(1));sv2=TModel.Vertices(slave.Vertices(2));
                    if(abs(v1.Y-sv1.Y)<=error && abs(v1.X-sv1.X)<=error && abs(v2.Y-sv2.Y)<=error && abs(v2.X-sv2.X)<=error),pair_index=slave.Index;
                    elseif(abs(v1.Y-sv2.Y)<=error && abs(v1.X-sv2.X)<=error && abs(v2.Y-sv1.Y)<=error && abs(v2.X-sv1.X)<=error),pair_index=-slave.Index;
                    end
                end
    end
end
function [pair_index] = findFacetPeriodicPair(TModel,facet,boundary),slave_Boundary=TModel.Boundaries(boundary.Param);error=1000*eps;
    v1=TModel.Vertices(facet.Vertices(1));v2=TModel.Vertices(facet.Vertices(2));v3=TModel.Vertices(facet.Vertices(3));
    switch abs(boundary.Axis)
        case 1,for ii=1:numel(slave_Boundary.Facets),slave=TModel.Facets(slave_Boundary.Facets(ii));
                    sv1=TModel.Vertices(slave.Vertices(1));sv2=TModel.Vertices(slave.Vertices(2));sv3=TModel.Vertices(slave.Vertices(3));
                    if(abs(sv1.Y-v1.Y)<=error && abs(sv2.Y-v2.Y)<=error && abs(sv3.Y-v3.Y)<=error && abs(sv1.Z-v1.Z)<=error && abs(sv2.Z-v2.Z)<=error && abs(sv3.Z-v3.Z)<=error),pair_index=slave.Index;
                    elseif(abs(sv1.Y-v3.Y)<=error && abs(sv2.Y-v1.Y)<=error && abs(sv3.Y-v2.Y)<=error && abs(sv1.Z-v3.Z)<=error && abs(sv2.Z-v1.Z)<=error && abs(sv3.Z-v2.Z)<=error),pair_index=slave.Index;
                    elseif(abs(sv1.Y-v2.Y)<=error && abs(sv2.Y-v3.Y)<=error && abs(sv3.Y-v1.Y)<=error && abs(sv1.Z-v2.Z)<=error && abs(sv2.Z-v3.Z)<=error && abs(sv3.Z-v1.Z)<=error),pair_index=slave.Index;
                    elseif(abs(sv1.Y-v1.Y)<=error && abs(sv2.Y-v3.Y)<=error && abs(sv3.Y-v2.Y)<=error && abs(sv1.Z-v1.Z)<=error && abs(sv2.Z-v3.Z)<=error && abs(sv3.Z-v2.Z)<=error),pair_index=-slave.Index;
                    elseif(abs(sv1.Y-v3.Y)<=error && abs(sv2.Y-v2.Y)<=error && abs(sv3.Y-v1.Y)<=error && abs(sv1.Z-v3.Z)<=error && abs(sv2.Z-v2.Z)<=error && abs(sv3.Z-v1.Z)<=error),pair_index=-slave.Index;
                    elseif(abs(sv1.Y-v2.Y)<=error && abs(sv2.Y-v1.Y)<=error && abs(sv3.Y-v3.Y)<=error && abs(sv1.Z-v2.Z)<=error && abs(sv2.Z-v1.Z)<=error && abs(sv3.Z-v3.Z)<=error),pair_index=-slave.Index;
                    end
               end
        case 2,for ii=1:numel(slave_Boundary.Facets),slave=TModel.Facets(slave_Boundary.Facets(ii));
                    sv1=TModel.Vertices(slave.Vertices(1));sv2=TModel.Vertices(slave.Vertices(2));sv3=TModel.Vertices(slave.Vertices(3));
                    if(abs(sv1.X-v1.X)<=error && abs(sv2.X-v2.X)<=error && abs(sv3.X-v3.X)<=error && abs(sv1.Z-v1.Z)<=error && abs(sv2.Z-v2.Z)<=error && abs(sv3.Z-v3.Z)<=error),pair_index=slave.Index;
                    elseif(abs(sv1.X-v3.X)<=error && abs(sv2.X-v1.X)<=error && abs(sv3.X-v2.X)<=error && abs(sv1.Z-v3.Z)<=error && abs(sv2.Z-v1.Z)<=error && abs(sv3.Z-v2.Z)<=error),pair_index=slave.Index;
                    elseif(abs(sv1.X-v2.X)<=error && abs(sv2.X-v3.X)<=error && abs(sv3.X-v1.X)<=error && abs(sv1.Z-v2.Z)<=error && abs(sv2.Z-v3.Z)<=error && abs(sv3.Z-v1.Z)<=error),pair_index=slave.Index;
                    elseif(abs(sv1.X-v1.X)<=error && abs(sv2.X-v3.X)<=error && abs(sv3.X-v2.X)<=error && abs(sv1.Z-v1.Z)<=error && abs(sv2.Z-v3.Z)<=error && abs(sv3.Z-v2.Z)<=error),pair_index=-slave.Index;
                    elseif(abs(sv1.X-v3.X)<=error && abs(sv2.X-v2.X)<=error && abs(sv3.X-v1.X)<=error && abs(sv1.Z-v3.Z)<=error && abs(sv2.Z-v2.Z)<=error && abs(sv3.Z-v1.Z)<=error),pair_index=-slave.Index;
                    elseif(abs(sv1.X-v2.X)<=error && abs(sv2.X-v1.X)<=error && abs(sv3.X-v3.X)<=error && abs(sv1.Z-v2.Z)<=error && abs(sv2.Z-v1.Z)<=error && abs(sv3.Z-v3.Z)<=error),pair_index=-slave.Index;
                    end
               end
        case 3,for ii=1:numel(slave_Boundary.Facets),slave=TModel.Facets(slave_Boundary.Facets(ii));
                    sv1=TModel.Vertices(slave.Vertices(1));sv2=TModel.Vertices(slave.Vertices(2));sv3=TModel.Vertices(slave.Vertices(3));
                    if(abs(sv1.Y-v1.Y)<=error && abs(sv2.Y-v2.Y)<=error && abs(sv3.Y-v3.Y)<=error && abs(sv1.X-v1.X)<=error && abs(sv2.X-v2.X)<=error && abs(sv3.X-v3.X)<=error),pair_index=slave.Index;
                    elseif(abs(sv1.Y-v3.Y)<=error && abs(sv2.Y-v1.Y)<=error && abs(sv3.Y-v2.Y)<=error && abs(sv1.X-v3.X)<=error && abs(sv2.X-v1.X)<=error && abs(sv3.X-v2.X)<=error),pair_index=slave.Index;
                    elseif(abs(sv1.Y-v2.Y)<=error && abs(sv2.Y-v3.Y)<=error && abs(sv3.Y-v1.Y)<=error && abs(sv1.X-v2.X)<=error && abs(sv2.X-v3.X)<=error && abs(sv3.X-v1.X)<=error),pair_index=slave.Index;
                    elseif(abs(sv1.Y-v1.Y)<=error && abs(sv2.Y-v3.Y)<=error && abs(sv3.Y-v2.Y)<=error && abs(sv1.X-v1.X)<=error && abs(sv2.X-v3.X)<=error && abs(sv3.X-v2.X)<=error),pair_index=-slave.Index;
                    elseif(abs(sv1.Y-v3.Y)<=error && abs(sv2.Y-v2.Y)<=error && abs(sv3.Y-v1.Y)<=error && abs(sv1.X-v3.X)<=error && abs(sv2.X-v2.X)<=error && abs(sv3.X-v1.X)<=error),pair_index=-slave.Index;
                    elseif(abs(sv1.Y-v2.Y)<=error && abs(sv2.Y-v1.Y)<=error && abs(sv3.Y-v3.Y)<=error && abs(sv1.X-v2.X)<=error && abs(sv2.X-v1.X)<=error && abs(sv3.X-v3.X)<=error),pair_index=-slave.Index;
                    end
               end
    end
end