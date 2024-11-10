function [] = UpdateCurrentBoundary(app),cla(app.UIAxes,'reset');
    boundary=app.TModel.Boundaries(app.CurrentBoundaryIndex);TModel=app.TModel;PlotBoundary_GUI(app,app.CurrentBoundaryIndex);
    if(boundary.Type=="PEC"),app.BoundaryTypeSelection.Value="Perfect Electric Conductor";app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;app.BoundaryTensorButton.Enable=false;
                             app.BoundaryTensorButton.Visible=false;app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
    elseif(boundary.Type=="PMC"),app.BoundaryTypeSelection.Value="Perfect Magnetic Conductor";app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;app.BoundaryTensorButton.Enable=false;
                                 app.BoundaryTensorButton.Visible=false;app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
    elseif(boundary.Type=="ABZ"),app.BoundaryTypeSelection.Value="Absorbing Boundary Condition";app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Enable=true;
        if(app.TModel.Frequency.NF>1),app.BoundaryDispersiveCheck.Enable=true;app.BoundaryDispersiveCheck.Visible=true;
            if(boundary.Dispersive==true),app.BoundaryDispersiveCheck.Value=true;end
        else,app.BoundaryDispersiveCheck.Enable=false;app.BoundaryDispersiveCheck.Visible=false;
        end,app.BoundaryTensorButton.Text="Ζ";app.PropagationConstantButton.Value=false;app.WaveImpedanceButton.Value=true;
    elseif(boundary.Type=="ABB"),app.BoundaryTypeSelection.Value="Absorbing Boundary Condition";app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Enable=true;
        if(app.TModel.Frequency.NF>1),app.BoundaryDispersiveCheck.Enable=true;app.BoundaryDispersiveCheck.Visible=true;
            if(boundary.Dispersive==true),app.BoundaryDispersiveCheck.Value=true;end
        else,app.BoundaryDispersiveCheck.Enable=false;app.BoundaryDispersiveCheck.Visible=false;
        end,app.BoundaryTensorButton.Text="β";app.PropagationConstantButton.Value=true;app.WaveImpedanceButton.Value=false;
    elseif(boundary.Type=="GRA"),app.BoundaryTypeSelection.Value="Graphene";app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Text="Conductivity";
        if(app.TModel.Frequency.NF>1),app.BoundaryDispersiveCheck.Enable=true;app.BoundaryDispersiveCheck.Visible=true;
            if(boundary.Dispersive==true),app.BoundaryDispersiveCheck.Value=true;end
        else,app.BoundaryDispersiveCheck.Enable=false;app.BoundaryDispersiveCheck.Visible=false;
        end
    elseif(boundary.Type=="IBC"),app.BoundaryTypeSelection.Value="Impedance Boundary Condition";app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Text="σ";
        if(app.TModel.Frequency.NF>1),app.BoundaryDispersiveCheck.Enable=true;app.BoundaryDispersiveCheck.Visible=true; app.BoundaryTensorButton.Visible=true;
            if(boundary.Dispersive==true),app.BoundaryDispersiveCheck.Value=true;end
        else,app.BoundaryDispersiveCheck.Enable=false;app.BoundaryDispersiveCheck.Visible=false;
        end
    elseif(boundary.Type=="PBC"),app.BoundaryTypeSelection.Value="Periodic Boundary Condition";app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;app.BoundaryTensorButton.Text="Periodic Pair";
                             app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Enable=true;app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                             if(~isempty(boundary.Param))

                             end
    elseif(boundary.Type=="CON"),app.BoundaryTypeSelection.Value="Continuity";app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;app.BoundaryTensorButton.Enable=false;
                             app.BoundaryTensorButton.Visible=false;app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                             if(~isempty(app.TModel.Boundaries(app.CurrentBoundaryIndex).ExcitationIndex)),app.BoundaryTensorButton.Enable=false;end
    elseif(boundary.Type=="POR"),app.BoundaryTypeSelection.Value="Port";
        SetPortBoundary(app);if(~isempty(app.TModel.Boundaries(app.CurrentBoundaryIndex).ExcitationIndex)),app.BoundaryTensorButton.Enable=false;end
    elseif(boundary.Type=="DIR"),app.BoundaryTypeSelection.Value="Dirichlet";
        SetPortBoundary(app);if(~isempty(app.TModel.Boundaries(app.CurrentBoundaryIndex).ExcitationIndex)),app.BoundaryTensorButton.Enable=false;end
    end
end



