function [] = UpdateCurrentBoundary_Excitation(app)
    boundary=app.TModel.Boundaries(app.CurrentBoundaryIndex);PlotBoundary_GUI(app,app.CurrentBoundaryIndex);
    switch boundary.Type
        case "PEC"
                app.BoundaryTypeSelection.Value="Perfect Electric Conductor";
                app.BoundaryTensorButton.Enable=false;app.BoundaryTensorButton.Visible=false;
                app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
        case "PMC"
                app.BoundaryTypeSelection.Value="Perfect Magnetic Conductor";
                app.BoundaryTensorButton.Enable=false;app.BoundaryTensorButton.Visible=false;
                app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
        case "CON"
                app.BoundaryTypeSelection.Value="Continuity";
                app.BoundaryTensorButton.Enable=false;app.BoundaryTensorButton.Visible=false;
                app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
        case "ABC"
                app.BoundaryTypeSelection.Value="Absorbing Boundary Condition";
                app.BoundaryTensorButton.Enable=false;app.BoundaryTensorButton.Visible=false;
                app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
        case "ABB"
                app.BoundaryTypeSelection.Value="Port ABC";
                app.BoundaryButtonGroup.Enable=true;app.BoundaryButtonGroup.Visible=true;
                app.WaveImpedanceButton.Text="Wave Impedance (Ζ)";app.PropagationConstantButton.Text="Propagation Constant (β)";
                app.WaveImpedanceButton.Value=false;app.PropagationConstantButton.Value=true;
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;
                if(app.TModel.Frequency.NF==1),app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                else,app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;
                end
                if(boundary.Dispersive),app.BoundaryDispersiveCheck.Value=true;else,app.BoundaryDispersiveCheck.Value=false;end
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Enable=true;
                app.BoundaryTensorButton.Text="β";
        case "ABZ"
                app.BoundaryTypeSelection.Value="Port ABC";
                app.BoundaryButtonGroup.Enable=true;app.BoundaryButtonGroup.Visible=true;
                app.WaveImpedanceButton.Text="Wave Impedance (Ζ)";app.PropagationConstantButton.Text="Propagation Constant (β)";
                if(app.TModel.Frequency.NF==1),app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                else,app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;
                end
                app.WaveImpedanceButton.Value=true;app.PropagationConstantButton.Value=false;
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;
                app.TensorCheckBox.Visible=true;app.TensorCheckBox.Enable=true;
                if(boundary.Dispersive),app.BoundaryDispersiveCheck.Value=true;else,app.BoundaryDispersiveCheck.Value=false;end
                if(boundary.Tensor),app.TensorCheckBox.Value=true;else,app.TensorCheckBox.Value=false;end
                app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Enable=true;
                app.BoundaryTensorButton.Text="Ζ";
        case "DIR"
                app.BoundaryTypeSelection.Value="Dirichlet";
                app.WaveImpedanceButton.Text="2D 1/2 EigenSolver";app.PropagationConstantButton.Text="Plane Wave";
                app.BoundaryButtonGroup.Visible=true;app.BoundaryButtonGroup.Enable=true;
                app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Enable=true;
                if(app.TModel.Frequency.NF==1),app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                else,app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;
                end
                if(boundary.Dispersive),app.BoundaryDispersiveCheck.Value=true;else,app.BoundaryDispersiveCheck.Value=false;end
                switch boundary.PortType
                    case 0
                        app.WaveImpedanceButton.Value=true;app.PropagationConstantButton.Value=false;
                        app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                        app.WaveImpedanceButton.Value=true;app.PropagationConstantButton.Value=false;
                        app.PortParamPanel.Enable=true;app.PortParamPanel.Visible=true;
                        switch(boundary.PortParamType)
                            case 0,app.PlaneWave_bButton.Value=true;app.PlaneWave_ZButton.Value=false;
                            case 1,app.PlaneWave_bButton.Value=false;app.PlaneWave_ZButton.Value=true;
                            case 2,app.PlaneWave_bButton.Value=false;app.PlaneWave_ZButton.Value=true;
                        end
                        app.BoundaryTensorButton.Text="Excitation";
                    case 1
                        app.WaveImpedanceButton.Value=false;app.PropagationConstantButton.Value=true;
                        app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                        app.WaveImpedanceButton.Value=true;app.PropagationConstantButton.Value=false;
                        app.PortParamPanel.Enable=true;app.PortParamPanel.Visible=true;
                        switch(boundary.PortParamType)
                            case 0,app.PlaneWave_bButton.Value=true;app.PlaneWave_ZButton.Value=false;app.BoundaryTensorButton.Text="β";
                                    
                            case 1,app.PlaneWave_bButton.Value=false;app.PlaneWave_ZButton.Value=true;app.TensorCheckBox.Enable=true;app.TensorCheckBox.Visible=true;app.TensorCheckBox.Value=false;
                                   app.BoundaryTensorButton.Text="Ζ";
                            case 2,app.PlaneWave_bButton.Value=false;app.PlaneWave_ZButton.Value=true;app.TensorCheckBox.Enable=true;app.TensorCheckBox.Visible=true;app.TensorCheckBox.Value=true;
                                   app.BoundaryTensorButton.Text="Ζ";
                        end
                        app.PlaneWavePanel.Visible=true;app.PlaneWavePanel.Enable=true;
                        set(app.PlaneWave_xField,"Value",boundary.PlaneWave(1));set(app.PlaneWave_yField,"Value",boundary.PlaneWave(2));set(app.PlaneWave_zField,"Value",boundary.PlaneWave(3));
                end
        case "POR"
                app.BoundaryTypeSelection.Value="Port";
                app.BoundaryButtonGroup.Visible=true;app.BoundaryButtonGroup.Enable=true;
                app.WaveImpedanceButton.Text="2D 1/2 EigenSolver";app.PropagationConstantButton.Text="Plane Wave";
                if(app.TModel.Frequency.NF==1),app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                else,app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;
                end
                if(boundary.Dispersive),app.BoundaryDispersiveCheck.Value=true;else,app.BoundaryDispersiveCheck.Value=false;end
                app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Enable=true;
                switch boundary.PortType
                    case 0
                        app.WaveImpedanceButton.Value=true;app.PropagationConstantButton.Value=false;
                        app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                        app.PortParamPanel.Enable=true;app.PortParamPanel.Visible=true;
                        switch(boundary.PortParamType)
                            case 0,app.PlaneWave_bButton.Value=true;app.PlaneWave_ZButton.Value=false;
                            case 1,app.PlaneWave_bButton.Value=false;app.PlaneWave_ZButton.Value=true;
                            case 2,app.PlaneWave_bButton.Value=false;app.PlaneWave_ZButton.Value=true;
                        end
                        app.BoundaryTensorButton.Text="Excitation";
                    case 1
                        app.WaveImpedanceButton.Value=false;app.PropagationConstantButton.Value=true;
                        app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                        app.PortParamPanel.Enable=true;app.PortParamPanel.Visible=true;
                        switch(boundary.PortParamType)
                            case 0,app.PlaneWave_bButton.Value=true;app.PlaneWave_ZButton.Value=false;app.BoundaryTensorButton.Text="β";
                                
                            case 1,app.PlaneWave_bButton.Value=false;app.PlaneWave_ZButton.Value=true;app.TensorCheckBox.Enable=true;app.TensorCheckBox.Visible=true;app.TensorCheckBox.Value=false;
                                   app.BoundaryTensorButton.Text="Ζ";
                            case 2,app.PlaneWave_bButton.Value=false;app.PlaneWave_ZButton.Value=true;app.TensorCheckBox.Enable=true;app.TensorCheckBox.Visible=true;app.TensorCheckBox.Value=true;
                                   app.BoundaryTensorButton.Text="Ζ";
                        end
                        app.PlaneWavePanel.Visible=true;app.PlaneWavePanel.Enable=true;
                        set(app.PlaneWave_xField,"Value",boundary.PlaneWave(1));set(app.PlaneWave_yField,"Value",boundary.PlaneWave(2));set(app.PlaneWave_zField,"Value",boundary.PlaneWave(3));
                end
        case "IBC"
                app.BoundaryTypeSelection.Value="Impedance Boundary Condition";
                app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                if(app.TModel.Frequency.NF==1),app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                else,app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;
                end
                if(boundary.Dispersive),app.BoundaryDispersiveCheck.Value=true;else,app.BoundaryDispersiveCheck.Value=false;end
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
                app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Text="Conductivity";
       case "PBC"
                app.BoundaryTypeSelection.Value="Absorbing Boundary Condition";
                app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
                app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Text="Periodic Pair";
        case "GRA"
                app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                if(app.TModel.Frequency.NF==1),app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                else,app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;
                end
                if(boundary.Dispersive),app.BoundaryDispersiveCheck.Value=true;else,app.BoundaryDispersiveCheck.Value=false;end
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
                app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Text="Conductivity";
    end
end
