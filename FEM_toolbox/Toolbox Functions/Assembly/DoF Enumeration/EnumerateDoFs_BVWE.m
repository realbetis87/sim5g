function [TModel,AssembledSystem] = EnumerateDoFs_BVWE(TModel,AssembledSystem,BoundaryIndices),UknownCounter=0;
    for kk=1:numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(kk));
        for ii=1:numel(boundary.Edges),edge=TModel.Edges(boundary.Edges(ii));
            if(isempty(edge.IndexE))
                if(isempty(edge.OnLine)),UknownCounter=UknownCounter+1;edge.IndexE=UknownCounter;
                elseif(edge.OnLine==0),UknownCounter=UknownCounter+1;edge.IndexE=UknownCounter;
                else,line=TModel.LineBoundaries(edge.OnLine);
                      if(line.Type~="PEC"),UknownCounter=UknownCounter+1;edge.IndexE=UknownCounter;else,edge.IndexE=0;end
                 end,TModel.Edges(boundary.Edges(ii))=edge;
            end
        end
    end,AssembledSystem.DimEt=UknownCounter;
    for kk=1:numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(kk));
            for ii=1:numel(boundary.Vertices),vertex=TModel.Vertices(boundary.Vertices(ii));
                if(isempty(vertex.IndexE))
                    if(isempty(vertex.OnLine)),UknownCounter=UknownCounter+1;vertex.IndexE=UknownCounter;
                    elseif(isscalar(vertex.OnLine))
                        if(vertex.OnLine==0),UknownCounter=UknownCounter+1;vertex.IndexE=UknownCounter;
                        else,line=TModel.LineBoundaries(vertex.OnLine);
                            if(line.Type~="PEC"),UknownCounter=UknownCounter+1;vertex.IndexE=UknownCounter;else,vertex.IndexE=0;end
                        end
                    else,pc=false;
                        for jj=1:numel(vertex.OnLine),line=TModel.LineBoundaries(vertex.OnLine(jj));if(line.Type=="PEC"),pc=true;end,end
                        if(pc),vertex.IndexE=0;else,UknownCounter=UknownCounter+1;vertex.IndexE=UknownCounter;end
                    end
                end
                TModel.Vertices(boundary.Vertices(ii))=vertex;
            end
    end,AssembledSystem.DimEn=UknownCounter;
end
