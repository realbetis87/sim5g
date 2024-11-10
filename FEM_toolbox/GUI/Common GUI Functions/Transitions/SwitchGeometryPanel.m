function [] = SwitchGeometryPanel(app)
        app.PreprocessingGrid.ColumnWidth={'1x','0x','0x','0x','0x'};
        app.CentralGrid.RowHeight={'1x','0x','0x','10x','1x'};
        app.ColumnMainGrid.ColumnWidth={'1x','4x','0x'};
        cla(app.UIAxes,'reset');PlotGeometry_GUI(app,0);
end

