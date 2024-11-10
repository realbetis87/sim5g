function [] = FrequencyDone_GUI(app),app.Step=2;
    app.PreprocessingGrid.ColumnWidth={'0x' '0x' '1x' '0x' '0x'};
    NewMessage(app,"Frequency Entered");
    if(app.TModel.Frequency.NF==1),app.DispersiveCheckBox.Enable=false;
    else,app.DispersiveCheckBox.Enable=true;
    end
    %---------------------- Disable Frequency Panel -----------------------
    app.StartStopPanel.Enable=false;
    app.SingleFrequencyCheckBox.Enable=false;app.SingleFrequencyField.Enable=false;app.SingleFrequencyUnit.Enable=false;
    app.MFButtonGroup.Enable=false;app.VectorField.Enable=false;app.FrequencyDoneButton.Enable=false;app.FrequencyDoneButton.Visible=false;
    %----------------------------------------------------------------------
    app.FrequencyButton.FontColor=[0.47,0.67,0.19];app.DomainsButton.Enable=true;app.TModel.Cond=3;
     cla(app.UIAxes,'reset');UpdateCurrentDomain(app);
end


% 0 Geometry,1 Frequency,2 Domains,3 Boundaries,4 Assembly,5 Solution.