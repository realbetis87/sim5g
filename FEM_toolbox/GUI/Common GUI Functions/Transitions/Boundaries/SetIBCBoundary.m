function [] = SetIBCBoundary(app),app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Text="Zs";app.BoundaryButtonGroup.Visible=false;app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
         if(app.TModel.Frequency.NF>1),app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;else,app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;end
end

