function [] = GeometryLoaded_EigenModeModule(app)
    %-------------- Left Panel --------------------------------------------
    app.LoadSelection.Enable=false;app.LoadSelection.Visible=false;app.LoadGeometryButton.Text="Initialize FEM Structures";
    app.LoadGrid.RowHeight={'0x','1x','2x'};app.ImportPanel.Title="";app.ImportPanel.BorderType="none";
    %-------------- Central Panel -----------------------------------------
    app.CentralGrid.RowHeight={'1x','0x','0x','10x','1x'};app.Image.Visible=false;app.Image.Enable=false;
    PlotGeometry_GUI(app,0);app.UIAxes.Visible=true;
end