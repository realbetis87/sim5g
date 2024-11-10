function [] = FrequencySelectionDone_GUI(app)
    switch(app.FrequencyMode)
           case 1,app.TModel=AddModelFrequency(app.TModel,app.SingleFrequencyField.Value,app.SingleFrequencyUnit.Value);
           case 2,Start=app.StartF.Value;Step=app.StepF.Value;Stop=app.StopF.Value;Unit=app.MFUnit.Value;
                  app.TModel=AddModelFrequency(app.TModel,Start,Step,Stop,Unit);
           case 3,app.TModel=AddModelFrequency(app.TModel,str2num(app.VectorField.Value),app.MFUnit.Value);
    end,FrequencyDone_GUI(app);pause(0.1);
end

