function [] = SwitchBoundaryPanel(app)
        app.PreprocessingGrid.ColumnWidth={'0x','0x','0x','1x','0x'};
        app.CentralGrid.RowHeight={'0x','1x','0x','10x','1x'};
        app.ColumnMainGrid.ColumnWidth={'1x','4x','1x'};app.RightPanelGrid.ColumnWidth={'1x','0x'};
         cla(app.UIAxes,'reset');UpdateCurrentBoundary(app)
end

