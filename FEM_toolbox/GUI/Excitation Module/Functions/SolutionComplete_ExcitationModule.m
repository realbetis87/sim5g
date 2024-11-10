function [] =SolutionComplete_ExcitationModule(app)
    app.ResultsPanel.Enable=true;app.ResultsPanel.Visible=true;
    if(app.TModel.Frequency.NF>1)
        app.ResFreqGrid.Visible=true;
        for ii=1:app.TModel.Frequency.NF,app.ResFreq.Items{end+1}=num2str(app.TModel.Frequency.UFrequency(ii));end,app.UnitLabel_2.Text=app.TModel.Frequency.Unit;
    end
end