function [] = GeometryAndImportsDone_EigenModeModule(app)
         app.PreprocessingGrid.ColumnWidth={'0x','1x','0x','0x','0x'};
         app.LoadGeometryButton.Enable=false;app.LoadGeometryButton.Visible=false;
         app.FrequencyButton.Enable=true;
end