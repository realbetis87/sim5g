function [] = PlotFields_2DModalSolver(app)
    Nhor=get(app.HNField,"Value");Nver=get(app.VNField,"Value");
    if(app.TModel.Frequency.NF==1)
        index=app.ResSelection.ValueIndex;
        switch app.Equation
            case 1,Plot2D_EH_TFLF(app.TModel,app.Excitation_Res,Nhor,Nver,app.EigenVectors(:,index));
            case 2,Plot2D_EH_TFTF(app.TModel,app.Excitation_Res,Nhor,Nver,app.EigenVectors(:,index));
            case 3,Plot2D_E_TFLF(app.TModel,app.Excitation_Res,Nhor,Nver,app.EigenVectors(:,index));
        end
    else
        index=app.ResSelection.ValueIndex;FreqIndex=app.FreqSelection.ValueIndex;EigenVector=app.EigenVectors{FreqIndex};
        switch app.Equation
            case 1,Plot2D_EH_TFLF(app.TModel,app.Excitation_Res,Nhor,Nver,EigenVector(:,index));
            case 2,Plot2D_EH_TFTF(app.TModel,app.Excitation_Res,Nhor,Nver,EigenVector(:,index));
            case 3,Plot2D_E_TFLF(app.TModel,app.Excitation_Res,Nhor,Nver,EigenVector(:,index));
        end
    end
end

