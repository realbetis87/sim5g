function [] = BoundariesDone_ExcitationModule(app)
    app.PreprocessingGrid.ColumnWidth={'0x','0x','0x','0x','1x'};app.RightPanelGrid.ColumnWidth={'0x','1x'};
    app.BoundariesButton.FontColor=[0.47,0.67,0.19];app.AssemblyButton.Enable=true;
    app.Step=4;app.TModel.Cond=5;
end

