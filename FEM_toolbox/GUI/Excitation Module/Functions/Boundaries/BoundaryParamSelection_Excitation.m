function [] = BoundaryParamSelection_Excitation(app)
    boundary=app.TModel.Boundaries(app.CurrentBoundaryIndex);
    switch boundary.Type
        case "POR"
            if(app.WaveImpedanceButton.Value)
                app.BoundaryTensorButton.Text="Excitation";
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;app.PortParamPanel.Visible=false;
                app.PlaneWave_bButton.Value=true;app.PlaneWave_ZButton.Value=false;
                boundary.PortType=0;boundary.PortParamType=0;
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=boundary;
                
            else
                app.BoundaryTensorButton.Text="β";
                app.PlaneWavePanel.Visible=true;app.PlaneWavePanel.Enable=true;
                set(app.PlaneWave_xField,"Value",0);set(app.PlaneWave_yField,"Value",0);set(app.PlaneWave_zField,"Value",0);
                app.PlaneWave_bButton.Value=true;app.PlaneWave_ZButton.Value=false;
                app.TModel.Boundaries(app.CurrentBoundaryIndex).PortType=1;
                app.TModel.Boundaries(app.CurrentBoundaryIndex).PortParamType=0;
                boundary.PortType=1;boundary.PortParamType=0;
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=boundary;
            end
        case "DIR"
            if(app.WaveImpedanceButton.Value)
                app.BoundaryTensorButton.Text="Excitation";
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;
                app.PlaneWave_bButton.Value=true;app.PlaneWave_ZButton.Value=false;
                app.TModel.Boundaries(app.CurrentBoundaryIndex).PortType=0;
                app.TModel.Boundaries(app.CurrentBoundaryIndex).PortParamType=0;
                boundary.PortType=0;boundary.PortParamType=0;
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=boundary;
            else
                app.BoundaryTensorButton.Text="β";
                app.PlaneWavePanel.Visible=true;app.PlaneWavePanel.Enable=true;
                set(app.PlaneWave_xField,"Value",0);set(app.PlaneWave_yField,"Value",0);set(app.PlaneWave_zField,"Value",0);
                app.PlaneWave_bButton.Value=true;app.PlaneWave_ZButton.Value=false;
                boundary.PortType=1;boundary.PortParamType=0;
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=boundary;
            end
        case "ABZ"
            if(app.PropagationConstantButton.Value)
                app.BoundaryTensorButton.Text="β";app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                if(app.BoundaryDispersiveCheck.Value),app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_B(app.TModel.Frequency);
                else,app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_B();
                end
            end
        case "ABB"
            if(app.WaveImpedanceButton.Value)
                app.BoundaryTensorButton.Text="Z";app.TensorCheckBox.Visible=true;app.TensorCheckBox.Enable=true;app.TensorCheckBox.Value=false;
                if(app.BoundaryDispersiveCheck.Value),app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_Z(app.TModel.Frequency);
                else,app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_Z();
                end
            end
     end
end

