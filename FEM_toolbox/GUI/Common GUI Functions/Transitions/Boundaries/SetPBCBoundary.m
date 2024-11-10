function [] = SetPBCBoundary(app)
    app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Text="Periodic Pair";app.BoundaryButtonGroup.Visible=false;app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
    app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
end

