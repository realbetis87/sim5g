function [] = Solve_ExcitationModule(app),Assembly=app.TModel.Assembled;
    if(app.TModel.Frequency.NF==1)
        if(Assembly.IsDir),app.TModel.Solution.ExcitationVector=app.TModel.Solution.ExcitationVector-app.TModel.Assembled.Matrix_B*app.TModel.Solution.KnownExcitation;end
        if(Assembly.IsPort),app.TModel.Solution.ExcitationVector=app.TModel.Solution.ExcitationVector+app.TModel.Assembled.TPV*app.TModel.Solution.UknownExcitation;end
        app.TModel.Solution.SolutionVector=app.TModel.Assembled.Matrix_A\app.TModel.Solution.ExcitationVector;
        NewMessage(app,"Single Frequency Excitation Solved");
    else
        if(app.SolutionIndex==0)
            for ff = 1: app.TModel.Frequency.NF
                if(Assembly.IsDir),app.TModel.Solution.ExcitationVector(:,ff)=app.TModel.Solution.ExcitationVector(:,ff)-app.TModel.Assembled.Matrix_B{ff}*app.TModel.Solution.KnownExcitation(:,ff);end
                if(Assembly.IsPort),app.TModel.Solution.ExcitationVector(:,ff)=app.TModel.Solution.ExcitationVector(:,ff)+app.TModel.Assembled.TPV{ff}*app.TModel.Solution.UknownExcitation(:,ff);end
                app.TModel.Solution.SolutionVector(:,ff)=app.TModel.Assembled.Matrix_A{ff}\app.TModel.Solution.ExcitationVector(:,ff);
            end
        else
            if(Assembly.IsDir),app.TModel.Solution.ExcitationVector(:,app.SolutionIndex)=app.TModel.Solution.ExcitationVector(:,app.SolutionIndex)-app.TModel.Assembled.Matrix_B{app.SolutionIndex}*app.TModel.Solution.KnownExcitation(:,app.SolutionIndex);end
            if(Assembly.IsPort),app.TModel.Solution.ExcitationVector(:,app.SolutionIndex)=app.TModel.Solution.ExcitationVector(:,app.SolutionIndex)+app.TModel.Assembled.TPV{app.SolutionIndex}*app.TModel.Solution.UknownExcitation(:,app.SolutionIndex);end
            app.TModel.Solution.SolutionVector(:,app.SolutionIndex)=app.TModel.Assembled.Matrix_A{app.SolutionIndex}\app.TModel.Solution.ExcitationVector(:,app.SolutionIndex);
        end
    end
    app.TModel.Cond=6;
end
