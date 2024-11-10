function [toolboxModel] = Prepare_PlaneWaveExcitation(toolboxModel),UVector=[];KVector=[];
    boundaries=toolboxModel.Boundaries;
    if(any([boundaries.PortType]==1))
        if(toolboxModel.Frequency.NF==1)
            for ib = 1 : numel(boundaries),boundary=toolboxModel.Boundaries(ib);
                if(boundary.Type=="POR")
                    if(boundary.PortType==1),if(isempty(UVector)),UVector=zeros(max([toolboxModel.Facets.UknownIndex]),1);end
                        for kk=1:numel(boundary.Edges)
                            edge=toolboxModel.Edges(boundary.Edges(kk));v1=toolboxModel.Vertices(edge.Vertices(1));v2=toolboxModel.Vertices(edge.Vertices(2));
                            UVector(edge.UknownIndex)=boundary.PlaneWave'*[v2.X-v1.X;v2.Y -v1.Y;v2.Z - v1.Z];
                        end
                    end                    
                elseif(boundary.Type=="DIR")
                    if(isempty(KVector)),KVector=zeros(max([toolboxModel.Facets.KnownIndex]),1);end
                     for kk=1:numel(boundary.Edges)
                            edge=toolboxModel.Edges(boundary.Edges(kk));v1=toolboxModel.Vertices(edge.Vertices(1));v2=toolboxModel.Vertices(edge.Vertices(2));
                            KVector(edge.KnownIndex)=boundary.PlaneWave'*[v2.X-v1.X;v2.Y -v1.Y;v2.Z - v1.Z];
                     end
                end
            end
        else
           for ib = 1 : numel(boundaries),boundary=toolboxModel.Boundaries(ib);
                if(boundary.Type=="POR")
                    if(boundary.PortType==1),if(isempty(UVector)),UVector=cell(max([toolboxModel.Facets.UknownIndex]),toolboxModel.Frequency.NF);end
                        for kk=1:numel(boundary.Edges)
                            edge=toolboxModel.Edges(boundary.Edges(kk));v1=toolboxModel.Vertices(edge.Vertices(1));v2=toolboxModel.Vertices(edge.Vertices(2));
                            UVector(edge.UknownIndex)=boundary.PlaneWave'*[v2.X-v1.X;v2.Y -v1.Y;v2.Z - v1.Z];
                        end
                    end                    
                elseif(boundary.Type=="DIR")
                    if(isempty(KVector)),KVector=cell(max([toolboxModel.Facets.KnownIndex]),toolboxModel.Frequency.NF);end
                     for kk=1:numel(boundary.Edges)
                            edge=toolboxModel.Edges(boundary.Edges(kk));v1=toolboxModel.Vertices(edge.Vertices(1));v2=toolboxModel.Vertices(edge.Vertices(2));
                            KVector(edge.KnownIndex)=boundary.PlaneWave'*[v2.X-v1.X;v2.Y -v1.Y;v2.Z - v1.Z];
                     end
                end
            end
        end
    end
end

