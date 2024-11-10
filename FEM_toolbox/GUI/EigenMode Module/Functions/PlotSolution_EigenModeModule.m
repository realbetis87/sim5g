function [] = PlotSolution_EigenModeModule(app)
    if(app.TModel.Frequency.NF==1)
        eig_index=app.ResSelection.ValueIndex;X=app.TModel.Solution.EigenVectors(:,eig_index);eigen_value=app.TModel.Solution.EigenValues(eig_index);
        Nv=get(app.HNField,"Value");Nh=get(app.VNField,"Value");Pos=get(app.PositionField,"Value");
        AxisIndex=app.PlotAxisSelection.ValueIndex;ax=app.TModel.Assembled.PropagationAxis;
        switch AxisIndex
            case 1,position=[Pos;nan;nan;];
            case 2,position=[nan;Pos;nan;];
            case 3,position=[nan;nan;Pos;];
        end
    %PlotSlice_EigenModeEB(app.TModel,X,position,[Nv,Nh],ax,eigen_value);
        PlotPlane_EB_EigenMode(app.TModel,X,position,[Nv,Nh],ax,eigen_value);
    else
        eig_index=app.ResSelection.ValueIndex;eigenVectors=app.TModel.Solution.EigenVectors;eigenValues=app.TModel.Solution.EigenValues;
        freq_index=app.ResFreq.ValueIndex;eigenVectors=eigenVectors{freq_index};eigenValues=eigenValues{freq_index};
        X=eigenVectors(:,eig_index);eigen_value=eigenValues(eig_index);
        Nv=get(app.HNField,"Value");Nh=get(app.VNField,"Value");Pos=get(app.PositionField,"Value");
        AxisIndex=app.PlotAxisSelection.ValueIndex;ax=app.TModel.Assembled.PropagationAxis;
        switch AxisIndex
            case 1,position=[Pos;nan;nan;];
            case 2,position=[nan;Pos;nan;];
            case 3,position=[nan;nan;Pos;];
        end
    %PlotSlice_EigenModeEB(app.TModel,X,position,[Nv,Nh],ax,eigen_value);
        PlotPlane_EB_EigenMode(app.TModel,X,position,[Nv,Nh],ax,eigen_value);

    end
end

