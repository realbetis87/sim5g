function [] = PortBoundarySelected(app)
    if(app.WaveImpedanceButton.Value)
        EH_EigenModeSolverForm(app.TModel,app.CurrentBoundaryIndex);
    else
    end
end

