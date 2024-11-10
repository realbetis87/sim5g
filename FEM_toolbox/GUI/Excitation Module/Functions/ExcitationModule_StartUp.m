function [] = ExcitationModule_StartUp(app)
   app.ColumnMainGrid.ColumnWidth={'1x','4x','0x'};
   %-------------------- Left Panel ---------------------------------------
   app.PreprocessingGrid.ColumnWidth={'1x','0x','0x','0x','0x'};
   app.RightPanelGrid.ColumnWidth={'1x','0x'};
   app.SolutionPanel.Enable=false;app.SolutionPanel.Visible=false;
   app.ResultsPanel.Enable=false;app.ResultsPanel.Visible=false;
   
end

