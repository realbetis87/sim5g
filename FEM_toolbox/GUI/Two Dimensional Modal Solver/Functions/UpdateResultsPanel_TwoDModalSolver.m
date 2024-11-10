function [] = UpdateResultsPanel_TwoDModalSolver(app)
    app.ResultsPanel.Enable=true;app.ResultsPanel.Visible=true;
    if(app.TModel.Frequency.NF==1),app.FreqSelection.Visible=false;app.FreqSelection.Enable=false;app.FreqUnit.Visible=false;app.ResSelection.Items={};
        for ii=1:numel(app.EigenValues),app.ResSelection.Items{end+1}=num2str(app.EigenValues(ii));end
           
    else,app.FreqSelection.Visible=true;app.FreqSelection.Enable=true;app.FreqUnit.Visible=true;
           for ii=1:app.TModel.Frequency.NF,app.FreqSelection.Items{end+1}=num2str(app.TModel.Frequency.UFrequency(ii));end,app.FreqUnit.Text=app.TModel.Frequency.Unit;
            value = app.FreqSelection.ValueIndex;eigenValues=app.EigenValues{value};app.ResSelection.Items={};
            for ii=1:numel(eigenValues)
                app.ResSelection.Items{end+1}=num2str(eigenValues(ii));
            end
    end
end

