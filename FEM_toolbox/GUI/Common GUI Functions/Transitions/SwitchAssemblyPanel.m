function [] = SwitchAssemblyPanel(app)
        app.PreprocessingGrid.ColumnWidth={'0x','0x','0x','0x','1x'};
        app.CentralGrid.RowHeight={'0x','1x','0x','10x','1x'};
        app.ColumnMainGrid.ColumnWidth={'1x','4x','1x'};app.RightPanelGrid.ColumnWidth={'0x','1x'};
         cla(app.UIAxes,'reset');PlotGeometry_GUI(app,0);
end

