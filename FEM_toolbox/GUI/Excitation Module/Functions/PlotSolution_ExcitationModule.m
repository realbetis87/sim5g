function [] = PlotSolution_ExcitationModule(app)
    if(app.TModel.Frequency.NF==1)
        X=app.TModel.Solution.SolutionVector;
        Nv=get(app.HNField,"Value");Nh=get(app.VNField,"Value");Pos=get(app.PositionField,"Value");
        AxisIndex=app.PlotAxisSelection.ValueIndex;
        switch AxisIndex
            case 1,position=[Pos;nan;nan;];
            case 2,position=[nan;Pos;nan;];
            case 3,position=[nan;nan;Pos;];
        end
        PlotPlane_EB(app.TModel,X,position,[Nv,Nh]);
    else
        freq_index=app.ResFreq.ValueIndex;X=app.TModel.Solution.SolutionVector(:,freq_index);
        Nv=get(app.HNField,"Value");Nh=get(app.VNField,"Value");Pos=get(app.PositionField,"Value");
        AxisIndex=app.PlotAxisSelection.ValueIndex;
        switch AxisIndex
            case 1,position=[Pos;nan;nan;];
            case 2,position=[nan;Pos;nan;];
            case 3,position=[nan;nan;Pos;];
        end
        PlotPlane_EB(app.TModel,X,position,[Nv,Nh]);
    end
end

