function [] = BoundaryDispersiveCheck_EigenModeModule(app)
            value = app.BoundaryDispersiveCheck.Value;boundary = app.BoundaryTypeSelection.Value;
            if(value)
                if(boundary=="Impedance Boundary Condition"),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).IBC(app.TModel.Frequency);
                elseif(boundary=="Absorbing Bondary Condition - Impedance")
                    if(app.TModel.Boundaries(app.CurrentBoundaryIndex).Tensor),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_ZT(app.TModel.Frequency);
                    else,app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_Z(app.TModel.Frequency);
                    end
                 elseif(boundary=="Graphene"),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).Graphene(app.TModel.Frequency);
                end
            else
                if(boundary=="Impedance Boundary Condition"),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).IBC();
                elseif(boundary=="Absorbing Bondary Condition - Impedance")
                    if(app.TModel.Boundaries(app.CurrentBoundaryIndex).Tensor),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_ZT();
                    else,app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).ABC_Z();
                    end
                elseif(boundary=="Graphene"),app.TModel.Boundaries(app.CurrentBoundaryIndex)=app.TModel.Boundaries(app.CurrentBoundaryIndex).Graphene();
                end
            end   
end