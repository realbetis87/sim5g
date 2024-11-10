function [TModel,AssembledSystem] = EnumerateDoFs_TFTF(TModel,AssembledSystem,BoundaryIndices),UknownCounter=0;
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
        for ii=1:numel(boundary.Edges),edge=TModel.Edges(boundary.Edges(ii));
            if(isempty(edge.IndexH))
                if(isempty(edge.OnLine)),UknownCounter=UknownCounter+1;edge.IndexH=UknownCounter;
                elseif(edge.OnLine==0),UknownCounter=UknownCounter+1;edge.IndexH=UknownCounter;
                else,line=TModel.LineBoundaries(edge.OnLine);
                      if(line.Type~="PMC"),UknownCounter=UknownCounter+1;edge.IndexH=UknownCounter;else,edge.IndexH=0;end
                 end,TModel.Edges(boundary.Edges(ii))=edge;
            end
        end
    end,AssembledSystem.DimHt=UknownCounter;
end
