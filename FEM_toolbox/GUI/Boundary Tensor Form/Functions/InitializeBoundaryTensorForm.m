function [] = InitializeBoundaryTensorForm(app),
    switch app.type
        case 1,app.boundary=app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex);
               app.BoundaryLabel.Text="Asorbing Boundary Condition";app.BoundaryTensorLabel.Text="Propagation Constant";app.ConductivityLabel.Text="β";app.boundary=app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex);
               if(app.boundary.Dispersive),app.Grid.RowHeight={'1x','1x','1x','0x','1x','0x','3x','1x'};set(app.TensorField,"Value",num2str(app.boundary.Param(1)));app.mode=1;
                for ii=1:app.mainApp.TModel.Frequency.NF,app.DropDown.Items{ii}=num2str(app.mainApp.TModel.Frequency.UFrequency(ii));end,app.UnitLabel.Text=app.mainApp.TModel.Frequency.Unit;app.DropDown.ValueIndex=1;
               else,app.Grid.RowHeight={'1x','1x','0x','0x','1x','0x','4x','1x'};set(app.TensorField,"Value",num2str(app.boundary.Param));app.mode=0;
               end
        case 2,boundary=app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex);app.boundary=app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex);
               if(boundary.Dispersive && boundary.Tensor),app.Grid.RowHeight={'1x','1x','1x','0x','0x','4x','0x','1x'};
                    for ii=1:app.mainApp.TModel.Frequency.NF,app.DropDown.Items{ii}=num2str(app.mainApp.TModel.Frequency.UFrequency(ii));end,app.UnitLabel.Text=app.mainApp.TModel.Frequency.Unit;app.DropDown.ValueIndex=1;
                    app.BoundaryLabel.Text="Absorbing Boundary Condition";
                    app.BoundaryTensorLabel.Text="Wave Impedace Z";app.mode=3;tensor=app.boundary.Param{1};
                    set(app.XXField,"Value",num2str(tensor(1,1)));set(app.XYField,"Value",num2str(tensor(1,2)));set(app.XZField,"Value",num2str(tensor(1,3)));
		            set(app.YXField,"Value",num2str(tensor(2,1)));set(app.YYField,"Value",num2str(tensor(2,2)));set(app.YZField,"Value",num2str(tensor(2,3)));
		            set(app.ZXField,"Value",num2str(tensor(3,1)));set(app.ZYField,"Value",num2str(tensor(3,2)));set(app.ZZField,"Value",num2str(tensor(3,3)));
                elseif(boundary.Dispersive)
                    app.Grid.RowHeight={'1x','1x','1x','0x','1x','0x','3x','1x'};
                    app.BoundaryLabel.Text="Absorbing Boundary Condition";
                    app.BoundaryTensorLabel.Text="Wave Impedace Z";
                    app.ConductivityLabel.Text="Z(Ω)";app.mode=1;
                    set(app.TensorField,"Value",num2str(app.boundary.Param(1)));
                    for ii=1:app.mainApp.TModel.Frequency.NF,app.DropDown.Items{ii}=num2str(app.mainApp.TModel.Frequency.UFrequency(ii));end,app.UnitLabel.Text=app.mainApp.TModel.Frequency.Unit;app.DropDown.ValueIndex=1;
                elseif(boundary.Tensor)
                    app.Grid.RowHeight={'1x','1x','0x','0x','0x','4x','1x','1x'};
                    app.BoundaryLabel.Text="Absorbing Boundary Condition";
                    app.BoundaryTensorLabel.Text="Wave Impedace Z";
                    app.ConductivityLabel.Text="Z(Ω)";app.mode=2;
		            tensor=app.boundary.Param;
                    set(app.XXField,"Value",num2str(tensor(1,1)));set(app.XYField,"Value",num2str(tensor(1,2)));set(app.XZField,"Value",num2str(tensor(1,3)));
		            set(app.YXField,"Value",num2str(tensor(2,1)));set(app.YYField,"Value",num2str(tensor(2,2)));set(app.YZField,"Value",num2str(tensor(2,3)));
		            set(app.ZXField,"Value",num2str(tensor(3,1)));set(app.ZYField,"Value",num2str(tensor(3,2)));set(app.ZZField,"Value",num2str(tensor(3,3)));
                else,app.Grid.RowHeight={'1x','1x','0x','0x','1x','0x','4x','1x'};
                     app.BoundaryLabel.Text="Absorbing Boundary Condition";
                     app.BoundaryTensorLabel.Text="Wave Impedace Z";
                     app.ConductivityLabel.Text="Z(Ω)";app.mode=0;
		             set(app.TensorField,"Value",num2str(app.boundary.Param));
                end
        case 3,app.boundary=app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex);
               app.BoundaryLabel.Text="Graphene Boundary Condition";app.boundary=app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex);
               app.BoundaryTensorLabel.Text="Graphene Conductivity";app.ConductivityLabel.Text="σ_{G}";
               if(app.boundary.Dispersive),app.Grid.RowHeight={'1x','1x','1x','0x','1x','0x','3x','1x'};set(app.TensorField,"Value",num2str(app.boundary.Param(1)));app.mode=1;
                     for ii=1:app.mainApp.TModel.Frequency.NF,app.DropDown.Items{ii}=num2str(app.mainApp.TModel.Frequency.UFrequency(ii));end,app.UnitLabel.Text=app.mainApp.TModel.Frequency.Unit;app.DropDown.ValueIndex=1;
               else,app.Grid.RowHeight={'1x','1x','0x','0x','1x','0x','4x','1x'};set(app.TensorField,"Value",num2str(app.boundary.Param));app.mode=0;
               end
        case 4,app.boundary=app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex);
               app.BoundaryLabel.Text="Impedance Boundary Condition";app.BoundaryTensorLabel.Text="Complex Surface Impedance";app.ConductivityLabel.Text="Z_{s}";app.boundary=app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex);
               if(~app.boundary.Dispersive),app.Grid.RowHeight={'1x','1x','0x','0x','1x','0x','4x','1x'};set(app.TensorField,"Value",num2str(app.boundary.Param));app.mode=0;
               else,app.Grid.RowHeight={'1x','1x','1x','0x','1x','0x','3x','1x'};set(app.TensorField,"Value",num2str(app.boundary.Param(1)));app.mode=1;
                   for ii=1:app.mainApp.TModel.Frequency.NF,app.DropDown.Items{ii}=num2str(app.mainApp.TModel.Frequency.UFrequency(ii));end,app.UnitLabel.Text=app.mainApp.TModel.Frequency.Unit;app.DropDown.ValueIndex=1;
               end
        case 5,app.boundary=app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex);
               app.Grid.RowHeight={'1x','0x','0x','3x','0x','0x','4x','1x'};app.boundary=app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex);app.mode=4;
               app.BoundaryLabel.Text="Periodic Boundary Condition";boundary=app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex);
               if(boundary.Master),app.MasterBoundaryLabel_1.Text="Boundary " + num2str(boundary.Index);app.MasterLabel.Text="Master";
                                   app.MasterBoundaryLabel_2.Text="Boundary " + num2str(boundary.Index);
                                   if(isempty(boundary.Param)),app.PBCGrid.RowHeight={'1x','1x','0x'};
                                         slaves=Check4Slave(app,boundary);for ii=1:numel(slaves),app.SlaveSelection.Items{ii}= num2str(slaves(ii));end
                                   else,app.SlaveBoundaryLabel.Text=num2str(boundary.Param);app.PBCGrid.RowHeight={'1x','0x','1x'};
                                   end
                else,app.MasterBoundaryLabel_1.Text="Boundary " + num2str(boundary.Index);app.MasterLabel.Text="Slave";app.SlaveLabel.Text="Master";app.PBCGrid.RowHeight={'1x','0x','1x'};app.SlaveBoundaryLabel.Text=num2str(boundary.Param);
                end
        case 6,app.BoundaryLabel.Text="Graphene Boundary Condition";app.boundary=app.mainApp.TModel.LineBoundaries(app.mainApp.selectedLine);
               app.BoundaryTensorLabel.Text="Graphene Conductivity";app.ConductivityLabel.Text="σ_{G}";
               if(app.boundary.Dispersive),app.Grid.RowHeight={'1x','1x','1x','0x','1x','0x','3x','1x'};set(app.TensorField,"Value",num2str(app.boundary.Param(1)));app.mode=1;
                   for ii=1:app.mainApp.TModel.Frequency.NF,app.DropDown.Items{ii}=num2str(app.mainApp.TModel.Frequency.UFrequency(ii));end,app.UnitLabel.Text=app.mainApp.TModel.Frequency.Unit;app.DropDown.ValueIndex=1;
               else,app.Grid.RowHeight={'1x','1x','1x','0x','1x','0x','3x','1x'};set(app.TensorField,"Value",num2str(app.boundary.Param));app.mode=0;
               end
     end
end

function [slave] = Check4Slave(app,boundary),slave=[];
    for ii=1:numel(app.mainApp.TModel.Boundaries),slave_can=app.mainApp.TModel.Boundaries(ii);
        if(slave_can.Axis==-boundary.Axis),slave(end+1)=ii;end
    end
end
