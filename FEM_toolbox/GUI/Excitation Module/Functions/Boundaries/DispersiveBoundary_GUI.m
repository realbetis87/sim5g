function [] = DispersiveBoundary_GUI(app)
    if(app.BoundaryDispersiveCheck.Value)
            app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).Disp(app.TModel.Frequency);
    else,   app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).NonDisp();
    end
end

