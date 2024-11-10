function [] = SetDIRBoundary(app)
    app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Text="Excitation";
    app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
    app.WaveImpedanceButton.Text="Import Field Distribution";app.PropagationConstantButton.Text="2D 1/2 EigenSolver";app.BoundaryButtonGroup.Enable=true;app.BoundaryButtonGroup.Visible=true;
    if(app.TModel.Frequency.NF>1),app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;else,app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;end
end

