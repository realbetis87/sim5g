function [] = ToolboxModelLoaded_EigenModeModule(app)
    switch app.TModel.Cond
        case 0,ToolboxGeometryLoaded(app);
        case 1,ToolboxFEMStructuresLoaded(app);
        case 2,ToolboxFrequencyLoaded(app);
        case 3,ToolboxDomainsLoaded(app);
        case 4,ToolboxBoundariesLoaded(app);
        case 5,ToolboxAssemblyLoaded(app);
        case 6,ToolboxSolutionLoaded(app);
    end
end

function [] = ToolboxGeometryLoaded(app)
    app.GeometryLoaded=true;GeometryLoaded_EigenModeModule(app)
end
function [] = ToolboxFEMStructuresLoaded(app),ToolboxGeometryLoaded(app);app.InitializationsComplete=true;InitializationsComplete_EigenModeModule(app);app.FrequencyButton.Enable=true;end
function [] = ToolboxFrequencyLoaded(app),ToolboxFEMStructuresLoaded(app);
    if(app.TModel.Frequency.NF==1)
        app.SingleFrequencyCheckBox.Value=true;app.SingleFrequencyCheckBox.Enable=false;app.SingleFrequencyCheckBox.Visible=true;
        app.MFPanel.Visible=false;app.MFPanel.Enable=false;app.StartStopPanel.Visible=false;app.StartStopPanel.Enable=false;app.VectorPanel.Visible=false;app.VectorPanel.Enable=false;
        app.FrequencyDoneButton.Visible=false;app.FrequencyDoneButton.Enable=false;
        app.FrequencyButton.FontColor=[0.47,0.67,0.19];app.DomainsButton.Enable=true;app.TModel.Cond=3;
        set(app.SingleFrequencyField,"Value",app.TModel.Frequency.UFrequency);
        switch app.TModel.Frequency.Unit
            case "Hz",set(app.SingleFrequencyUnit,"ValueIndex",1);
            case "KHz",set(app.SingleFrequencyUnit,"ValueIndex",2);  
            case "MHz",set(app.SingleFrequencyUnit,"ValueIndex",3);
            case "GHz",set(app.SingleFrequencyUnit,"ValueIndex",4);
            case "THz",set(app.SingleFrequencyUnit,"ValueIndex",5);
        end
        app.SingleFrequencyPanel.Enable=false;app.SingleFrequencyPanel.Visible=true;
    else
        app.SingleFrequencyCheckBox.Value=false;app.SingleFrequencyCheckBox.Enable=false;app.SingleFrequencyCheckBox.Visible=true;
        app.FrequencyButton.FontColor=[0.47,0.67,0.19];app.DomainsButton.Enable=true;app.TModel.Cond=3;
        app.SingleFrequencyPanel.Visible=false;app.SingleFrequencyPanel.Enable=false;
        app.FrequencyDoneButton.Visible=false;app.FrequencyDoneButton.Enable=false;
        app.MFPanel.Visible=true;app.MFPanel.Enable=false;app.MFButtonGroup.Enable=false;app.MFButtonGroup.Visible=true;
        if(~isempty(app.TModel.Frequency.Start))
            app.StartStopPanel.Visible=true;app.StartStopPanel.Enable=false;
            app.VectorPanel.Visible=false;app.VectorPanel.Enable=false;
            app.StartStepStopButton.Value=true;app.MFButtonGroup.Visible=true;app.MFButtonGroup.Enable=false;
            switch app.TModel.Frequency.Unit
                case "Hz",set(app.MFUnit,"ValueIndex",1);
                case "KHz",set(app.MFUnit,"ValueIndex",2);  
                case "MHz",set(app.MFUnit,"ValueIndex",3);
                case "GHz",set(app.MFUnit,"ValueIndex",4);
                case "THz",set(app.MFUnit,"ValueIndex",5);
            end
            set(app.StartF,"Value",app.TModel.Frequency.Start);set(app.StepF,"Value",app.TModel.Frequency.Increment);set(app.StopF,"Value",app.TModel.Frequency.Stop);
        else
            app.StartStopPanel.Visible=false;app.StartStopPanel.Enable=false;
            app.VectorPanel.Visible=true;app.VectorPanel.Enable=false;
            app.StartStepStopButton.Value=true;app.MFButtonGroup.Visible=true;app.MFButtonGroup.Enable=false;
            switch app.TModel.Frequency.Unit
                case "Hz",set(app.MFUnit,"ValueIndex",1);
                case "KHz",set(app.MFUnit,"ValueIndex",2);  
                case "MHz",set(app.MFUnit,"ValueIndex",3);
                case "GHz",set(app.MFUnit,"ValueIndex",4);
                case "THz",set(app.MFUnit,"ValueIndex",5);
            end
            string=num2str(app.TModel.Frequency.UFrequency(1))+";";for ii=2:app.TModel.Frequency.NF,string=string+num2str(app.TModel.Frequency.UFrequency(ii))+";";end
            set(app.VectorField,"Value",string);
        end
    end
end
function [] = ToolboxDomainsLoaded(app),ToolboxFrequencyLoaded(app);UpdateDomainsPanel(app);app.DomainsButton.FontColor=[0.47,0.67,0.19];app.BoundariesButton.Enable=true;UpdateCurrentDomain(app);
end
function [] = ToolboxBoundariesLoaded(app),ToolboxDomainsLoaded(app);UpdateBoundariesPanel(app);app.BoundariesButton.FontColor=[0.47,0.67,0.19];app.AssemblyButton.Enable=true;UpdateCurrentBoundary(app);
end
function [] =ToolboxAssemblyLoaded(app),ToolboxBoundariesLoaded(app);
    if(app.TModel.Assembled.EigenValue=="k"),app.KCheck.Value=true;app.NCheck.Value=false;
    else,app.KCheck.Value=false;app.NCheck.Value=true;
    end
    switch app.TModel.Assembled.PropagationAxis
        case "x",app.PropagationAxisSelection.ValueIndex=1;
        case "y",app.PropagationAxisSelection.ValueIndex=2;
        case "z",app.PropagationAxisSelection.ValueIndex=3;
    end
    app.AssemblyButton.FontColor=[0.47,0.67,0.19];
    app.ExportMatricesButton.Enable=true;app.ExportMatricesButton.Visible=true;
    app.MatrixSparsitiesButton.Enable=true;app.MatrixSparsitiesButton.Visible=true;
    app.NE_Info.Text=num2str(app.TModel.Assembled.NE);
    app.NB_Info.Text=num2str(app.TModel.Assembled.NB);
    app.N_Info.Text=num2str(app.TModel.Assembled.N);
    AssemblyComplete_EigenModeModule(app);pause(0.1);
end
function [] = ToolboxSolutionLoaded(app)
ToolboxAssemblyLoaded(app);
UpdateResultsPanel_EigenModeModule(app);
end