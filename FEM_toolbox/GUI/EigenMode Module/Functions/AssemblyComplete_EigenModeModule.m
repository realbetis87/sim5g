function [] = AssemblyComplete_EigenModeModule(app)
app.SolutionPanel.Enable=true;app.SolutionPanel.Visible=true;app.MatrixSparsitiesButton.Enable=true;
app.ExportMatricesButton.Enable=true;
    app.AssemblyButton.FontColor=[0.47,0.67,0.19];
    if(app.TModel.Frequency.NF==1)
        app.SolFreqGrid.Visible=false;
        if(isempty(app.TModel.Solution)),app.TModel.Solution=Solution();app.TModel.Solution.Type="EigenMode";end
    else,app.SolFreqGrid.Visible=true;app.SolutionFreqSelection.Items{1}='Solve All Frequencies';
        for ii=2:app.TModel.Frequency.NF+1,app.SolutionFreqSelection.Items{ii}=num2str(app.TModel.Frequency.UFrequency(ii-1));end
        app.UnitLabel.Text=app.TModel.Frequency.Unit;app.SolutionFreqSelection.ValueIndex=1;
        app.EigenValues=cell(app.TModel.Frequency.NF,1);app.EigenVectors=cell(app.TModel.Frequency.NF,1);
        if(isempty(app.TModel.Solution)),app.TModel.Solution=Solution();app.TModel.Solution.Type="EigenMode";end
        app.TModel.Solution.EigenVectors=cell(app.TModel.Frequency.NF,1);app.TModel.Solution.EigenValues=cell(app.TModel.Frequency.NF,1);
    end
end
