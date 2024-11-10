function [] = SelectEigenvalue_2DModalSolver(app),index=app.ResSelection.ValueIndex;
    if(app.TModel.Frequency.NF==1)
        app.Excitation_Res.Vector=app.EigenVectors(:,index);
        if(app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex).PortParamType==0),app.Excitation_Res.Beta=app.EigenValues(index);
            for ii=1:numel(app.BoundaryIndices),app.TModel.Boundaries(app.BoundaryIndices(ii)).Param=app.EigenValues(index);end
        else,app.Excitation_Res.Zita=app.EigenValues(index);%Calculate Zita
        end
        app.DoneButton.Enable=true;
    else
        freq_index=app.FreqSelection.ValueIndex;EigenVector=app.EigenVectors{freq_index};EigenValues=app.EigenValues{freq_index};
        app.Excitation_Res.Vector{freq_index}=EigenVector(:,index);
        if(app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex).PortParamType==0),app.Excitation_Res.Beta(freq_index)=EigenValues(index);
        else,app.Excitation_Res.Zita(freq_index)=app.EigenValues(index);%Calculate Zita
        end
        NewMessage(app," Eigen Pair for frequency : " + num2str(app.TModel.Frequency.UFrequency(freq_index)) + app.TModel.Frequency.Unit + " defined");
        app.eigenvalues_selected(freq_index)=1;
        if(all([app.eigenvalues_selected]==1)),app.DoneButton.Enable=true;end
    end
end
