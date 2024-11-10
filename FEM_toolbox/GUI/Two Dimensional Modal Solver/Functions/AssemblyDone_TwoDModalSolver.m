function [] = AssemblyDone_TwoDModalSolver(app),app.SolutionPanel.Visible=true;app.SolutionPanel.Enable=true;
    if(app.TModel.Frequency.NF>1)
        app.SolFreq.Enable=true;app.SolFreq.Visible=true;
        app.SolFreq.Items{1}='Solve All Frequencies';
        for ii=1:app.TModel.Frequency.NF,app.SolFreq.Items{end+1}=num2str(app.TModel.Frequency.UFrequency(ii));end
        app.FreqUnit_2.Enable=true;app.FreqUnit_2.Visible=true;app.FreqUnit_2.Text=app.TModel.Frequency.Unit;
    else
        app.SolFreq.Enable=false;app.SolFreq.Visible=false;
        app.FreqUnit_2.Enable=false;app.FreqUnit_2.Visible=false;
    end
   
end

