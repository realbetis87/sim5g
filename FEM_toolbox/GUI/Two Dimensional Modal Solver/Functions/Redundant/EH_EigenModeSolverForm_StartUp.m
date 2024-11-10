function [] = EH_EigenModeSolverForm_StartUp(app)
    app.TModel=Export2DModalBoundary(app.TModel,app.BoundaryIndex);
    for ii=1:numel(app.TModel.Boundaries(app.BoundaryIndex).Lines),app.BoundarySelection.Items{end+1}=char(int2str(ii));end,app.selectedLine=1;
    Plot2DDomain_GUI(app);EquationSelection(app);
end

function [] = EquationSelection(app)
end