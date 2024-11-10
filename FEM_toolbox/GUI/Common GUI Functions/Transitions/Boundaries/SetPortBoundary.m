function [] = SetPortBoundary(app)
    if(app.WaveImpedanceButton.Value),app.BoundaryTensorButton.Text="Value";
    else,app.BoundaryTensorButton.Text="Excitation";
    end
    app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;
    app.BoundaryButtonGroup.Visible=false;app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
    app.WaveImpedanceButton.Text="2D 1/2 EigenSolver";app.PropagationConstantButton.Text="Plane Wave";app.BoundaryButtonGroup.Enable=true;app.BoundaryButtonGroup.Visible=true;
    if(app.TModel.Frequency.NF>1),app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;else,app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;end
end

