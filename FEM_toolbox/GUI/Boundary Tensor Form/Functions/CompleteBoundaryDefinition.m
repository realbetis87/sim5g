function [] = CompleteBoundaryDefinition(app)
    switch app.type
        case 1,app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex)=app.boundary;
        case 2,app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex)=app.boundary;
        case 3,app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex)=app.boundary;
        case 4,app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex)=app.boundary;
        case 5,app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex)=app.boundary;
               if(~isempty(app.boundary.Param)),slave=app.mainApp.TModel.Boundaries(app.boundary.Param);slave=slave.PBC("s");slave.Param=app.boundary.Index;app.mainApp.TModel.Boundaries(slave.Index)=slave;
              app.mainApp.BoundariesInfo{slave.Index,3}='PBC';app.mainApp.BoundariesTable.Data=app.mainApp.BoundariesInfo;
               end
        case 6,app.mainApp.TModel.LineBoundaries(app.mainApp.selectedLine)=app.boundary;
    end
end

