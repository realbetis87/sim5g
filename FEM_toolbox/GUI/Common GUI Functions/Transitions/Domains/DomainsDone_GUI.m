function [] = DomainsDone_GUI(app)
    app.DomainsButton.FontColor=[0.47,0.67,0.19];app.BoundariesButton.Enable=true;
    app.PreprocessingGrid.ColumnWidth={'0x','0x','0x','1x','0x'};
    UpdateBoundariesPanel(app);app.Step=3;
    app.IsotropicCheckBox.Enable=false;app.DispersiveCheckBox.Enable=false;app.BianisotropicCheckBox.Enable=false;
    app.DomainDoneButton.Enable=false;app.DomainDoneButton.Visible=false;app.CurrentDomainDone.Enable=false;app.CurrentDomainDone.Visible=false;app.TModel.Cond=4;
    app.CentralGrid.RowHeight={'0x','1x','0x','10x','1x'};app.PlotBoundaryButton.Enable=true;app.PlotBoundaryButton.Visible=true;
    app.PlotBoundaryPanel.Visible=true;app.PlotBoundaryPanel.Enable=true;app.PlotBoundarySelection.Visible=true;app.PlotBoundarySelection.Enable=true;
    app.ClearPlotBoundaryButton.Visible=true;app.ClearPlotBoundaryButton.Enable=true;
    NewMessage(app,"All Domains Defined");
    cla(app.UIAxes,'reset');PlotBoundary_GUI(app,app.CurrentBoundaryIndex);
end

