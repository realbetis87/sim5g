function [] = BoundaryTypeSelection_ExcitationModule(app,value)
    if(value=="Perfect Electric Conductor")
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).PEC();
                app.BoundariesInfo{app.CurrentBoundaryIndex,3}='PEC';
                app.BoundaryTensorButton.Enable=false;app.BoundaryTensorButton.Visible=false;
                app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
     elseif(value=="Perfect Magnetic Conductor")
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).PMC();
                app.BoundariesInfo{app.CurrentBoundaryIndex,3}='PMC';
                app.BoundaryTensorButton.Enable=false;app.BoundaryTensorButton.Visible=false;
                app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
     elseif(value=="Impedance Boundary Condition")
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).IBC();
                app.BoundariesInfo{app.CurrentBoundaryIndex,3}='IBC';
                app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;
                app.BoundaryTensorButton.Text="Z";app.BoundaryButtonGroup.Visible=false;
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                 if(app.TModel.Frequency.NF>1),app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;app.BoundaryDispersiveCheck.Value=false;
                 else,app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                 end
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
     elseif(value=="Absorbing Boundary Condition")
                app.BoundaryTensorButton.Enable=false;app.BoundaryTensorButton.Visible=false;
                app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC();
                app.BoundariesInfo{app.CurrentBoundaryIndex,3}='ABC';
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
     elseif(value=="Graphene")
                app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Text="Conductivity";
                app.BoundaryButtonGroup.Visible=false;app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                if(app.TModel.Frequency.NF>1),app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;app.BoundaryDispersiveCheck.Value=false;
                else,app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                end
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).Graphene();
                app.BoundariesInfo{app.CurrentBoundaryIndex,3}='GRA';
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
     elseif(value=="Continuity")
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).CON();
                app.BoundariesInfo{app.CurrentBoundaryIndex,3}='CON';
                app.BoundaryTensorButton.Enable=false;app.BoundaryTensorButton.Visible=false;
                app.BoundaryButtonGroup.Enable=false;app.BoundaryButtonGroup.Visible=false;
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
     elseif(value=="Periodic Boundary Condition")
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).PBC("m");
                app.BoundariesInfo{app.CurrentBoundaryIndex,3}='PBC';
                app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;app.BoundaryTensorButton.Text="Periodic Pair";
                app.BoundaryButtonGroup.Visible=false;app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
     elseif(value=="Port ABC"),SetABCBoundary(app);
                app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;
                app.BoundaryButtonGroup.Enable=true;app.BoundaryButtonGroup.Visible=true;
                app.WaveImpedanceButton.Text="Wave Impedance (Ζ)";app.PropagationConstantButton.Text="Propagation Constant (β)";
                 app.WaveImpedanceButton.Value=false;app.PropagationConstantButton.Value=true;
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_B();
                app.BoundariesInfo{app.CurrentBoundaryIndex,3}='ABB';app.BoundaryTensorButton.Text="β";
                app.TensorCheckBox.Visible=false;app.TensorCheckBox.Enable=false;
                if(app.TModel.Frequency.NF>1),app.BoundaryDispersiveCheck.Visible=true;app.BoundaryDispersiveCheck.Enable=true;app.BoundaryDispersiveCheck.Value=false;
                else,app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                end
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
     elseif(value=="Dirichlet")
                app.WaveImpedanceButton.Text="2D 1/2 EigenSolver";app.PropagationConstantButton.Text="Plane Wave";
                app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;
                app.BoundaryButtonGroup.Enable=true;app.BoundaryButtonGroup.Visible=true;
                app.WaveImpedanceButton.Value=true;app.PropagationConstantButton.Value=false;
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
                app.PlaneWave_bButton.Value=true;app.PlaneWave_ZButton.Value=false;
                app.BoundaryTensorButton.Text="Excitation";
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).Dirichlet();
                app.TModel.Boundaries(app.CurrentBoundaryIndex).PortType=0;
                app.TModel.Boundaries(app.CurrentBoundaryIndex).PortParamType=0;
                if(app.TModel.Frequency.NF>1),app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;app.BoundaryDispersiveCheck.Value=false;
                else,app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                end
                app.BoundariesInfo{app.CurrentBoundaryIndex,3}='DIR';
    elseif(value=="Port")
                app.WaveImpedanceButton.Text="2D 1/2 EigenSolver";app.PropagationConstantButton.Text="Plane Wave";
                app.BoundaryButtonGroup.Enable=true;app.BoundaryButtonGroup.Visible=true;
                app.WaveImpedanceButton.Value=true;app.PropagationConstantButton.Value=false;
                app.PlaneWavePanel.Visible=false;app.PlaneWavePanel.Enable=false;
                app.PortParamPanel.Visible=false;app.PortParamPanel.Enable=false;
                app.PlaneWave_bButton.Value=true;app.PlaneWave_ZButton.Value=false;
                app.BoundaryTensorButton.Enable=true;app.BoundaryTensorButton.Visible=true;
                app.BoundaryTensorButton.Text="Excitation";
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).Port();
                app.TModel.Boundaries(app.CurrentBoundaryIndex).PortType=0;
                app.TModel.Boundaries(app.CurrentBoundaryIndex).PortParamType=0;
                if(app.TModel.Frequency.NF>1),app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;app.BoundaryDispersiveCheck.Value=false;
                else,app.BoundaryDispersiveCheck.Visible=false;app.BoundaryDispersiveCheck.Enable=false;
                end
                app.BoundariesInfo{app.CurrentBoundaryIndex,3}='Port';
     end,app.BoundariesTable.Data=app.BoundariesInfo;
     pause(0.1);
end

