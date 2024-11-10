function [] = SetABCBoundary(app)
        app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;app.BoundaryButtonGroup.Enable=true;app.BoundaryButtonGroup.Visible=true;
        app.WaveImpedanceButton.Text="Wave Impedance (Ζ)";app.PropagationConstantButton.Text="Propagation Constant (β)";
        if(app.PropagationConstantButton.Value),app.BoundaryTensorButton.Text="β";app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;else,app.BoundaryTensorButton.Text="Z";app.TensorCheckBox.Visible=true;app.TensorCheckBox.Enable=true;end
        if(app.TModel.Frequency.NF>1),app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;else,app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;end
end

