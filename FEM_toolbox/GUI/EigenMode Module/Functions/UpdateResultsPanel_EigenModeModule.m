function [] = UpdateResultsPanel_EigenModeModule(app),app.ResultsPanel.Visible=true;app.ResultsPanel.Enable=true;
    if(app.TModel.Frequency.NF==1)
        app.ResultsGrid.RowHeight={'0x','1x','1x','1x','1x','1x','1x','1x'};app.ResFreq.Visible=false;app.ResFreq.Enable=false;
        app.ResSelection.Items={};
        for ii=1:numel(app.TModel.Solution.EigenValues),app.ResSelection.Items{end+1}=num2str(app.TModel.Solution.EigenValues(ii));end
    else
        app.ResultsGrid.RowHeight={'1x','1x','1x','1x','1x','1x','1x','1x'};app.ResSelection.Items={};
        for ii=1:app.TModel.Frequency.NF,app.ResFreq.Items{end+1}=num2str(app.TModel.Frequency.UFrequency(ii));end,app.UnitLabel_2.Text=app.TModel.Frequency.Unit;
        Index=app.ResFreq.ValueIndex;
        for ii=1:numel(app.TModel.Solution.EigenValues{Index}),eigenValues=app.TModel.Solution.EigenValues{Index};app.ResSelection.Items{end+1}=num2str(eigenValues(ii));end
    end

end

