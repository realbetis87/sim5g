function [] = BoundaryTypeSelection_GUI(app,value)
            if(value=="Perfect Electric Conductor"),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).PEC();
                app.BoundariesInfo{app.CurrentBoundaryIndex,3}='PEC';SetPECBoundary(app);
            elseif(value=="Perfect Magnetic Conductor"),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).PMC();
                app.BoundariesInfo{app.CurrentBoundaryIndex,3}='PMC';SetPMCBoundary(app);
            elseif(value=="Impedance Boundary Condition"),SetIBCBoundary(app);
                if(app.BoundaryDispersiveCheck.Value==false),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).IBC();
                else,app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).IBC(app.TModel.Frequency);
                end,app.BoundariesInfo{app.CurrentBoundaryIndex,3}='IBC';
            elseif(value=="Absorbing Boundary Condition"),SetABCBoundary_EigenModeModule(app);
                app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC();app.BoundariesInfo{app.CurrentBoundaryIndex,3}='ABC';
            elseif(value=="Graphene"),SetGRABoundary(app);
                if(app.BoundaryDispersiveCheck.Value==false),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).Graphene();
                else,app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).Graphene(app.TModel.Frequency);
                end,app.BoundariesInfo{app.CurrentBoundaryIndex,3}='GRA';
            elseif(value=="Continuity"),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).CON();SetPECBoundary(app);app.BoundariesInfo{app.CurrentBoundaryIndex,3}='CON';
            elseif(value=="Periodic Boundary Condition"),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).PBC("m");SetPBCBoundary(app);app.BoundariesInfo{app.CurrentBoundaryIndex,3}='PBC';
            elseif(value=="Port ABC"),SetABCBoundary(app);
                 if(app.PropagationConstantButton.Value)
                      if(app.BoundaryDispersiveCheck.Value==false),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_B();
                      else,app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_B(app.TModel.Frequency);
                      end,app.BoundariesInfo{app.CurrentBoundaryIndex,3}='ABB';
                 else
                     if(app.TensorCheckBox.Value)
                        if(app.BoundaryDispersiveCheck.Value==false),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_ZT();
                        else,app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_ZT(app.TModel.Frequency);
                        end,app.BoundariesInfo{app.CurrentBoundaryIndex,3}='ABZT';
                     else
                        if(app.BoundaryDispersiveCheck.Value==false),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_Z();
                        else,app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_Z(app.TModel.Frequency);
                        end,app.BoundariesInfo{app.CurrentBoundaryIndex,3}='ABZ';
                     end
                 end
            elseif(value=="Dirichlet"),SetPortBoundary(app);
                 if(app.BoundaryDispersiveCheck.Value==false),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).Dirichlet();
                 else,app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).Dirichlet(app.TModel.Frequency);
                end,app.BoundariesInfo{app.CurrentBoundaryIndex,3}='DIR';
                if(~isempty(app.TModel.Boundaries(app.CurrentBoundaryIndex).ExcitationIndex)),app.BoundaryTensorButton.Enable=false;end
            elseif(value=="Port"),SetPortBoundary(app);
                if(app.BoundaryDispersiveCheck.Value==false),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).Port();
                else,app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).Port(app.TModel.Frequency);
                end,app.BoundariesInfo{app.CurrentBoundaryIndex,3}='Port';
                if(~isempty(app.TModel.Boundaries(app.CurrentBoundaryIndex).ExcitationIndex)),app.BoundaryTensorButton.Enable=false;end
            end,app.BoundariesTable.Data=app.BoundariesInfo;pause(0.1);
end

