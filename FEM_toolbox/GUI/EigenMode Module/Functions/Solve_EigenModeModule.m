function [] = Solve_EigenModeModule(app),app.NumberOfEigenValues=get(app.EigenField,"Value");app.Shift=get(app.ShiftField,"Value");
    if(app.TModel.Frequency.NF>1)
         if(app.SolutionIndex==0)
                    tic;for ii=1:app.TModel.Frequency.NF,[V,D]=eigs(app.TModel.Assembled.Matrix_A{ii},app.TModel.Assembled.Matrix_B{ii},app.NumberOfEigenValues,app.Shift);
                            app.EigenVectors{ii}=V;app.EigenValues{ii}=diag(D);NewMessage(app,"Problem Solved for Frequency : " + num2str(app.TModel.Frequency.UFrequency)+" "+app.TModel.Frequency.Unit);
                            app.TModel.Solution.EigenVectors{ii}=V;app.TModel.Solution.EigenValues{ii}=diag(D);
                        end,SolveTime=toc;NewMessage(app,"Total Solution Time : " + num2str(SolveTime));pause(0.1);
         else
                    tic;[V,D]=eigs(app.TModel.Assembled.Matrix_A{app.SolutionIndex},app.TModel.Assembled.Matrix_B{app.SolutionIndex},app.NumberOfEigenValues,app.Shift);SolveTime=toc;
                    app.EigenVectors{app.SolutionIndex}=V;app.EigenValues{app.SolutionIndex}=diag(D);
                    app.TModel.Solution.EigenVectors{app.SolutionIndex}=V;app.TModel.Solution.EigenValues{app.SolutionIndex}=diag(D);
                    NewMessage(app,"Problem Solved for Frequency : " + num2str(app.TModel.Frequency.UFrequency)+" "+app.TModel.Frequency.Unit);
                    NewMessage(app,"Solution Time : " + num2str(SolveTime));pause(0.1);
         end
        
    else
        tic;[V,D]=eigs(app.TModel.Assembled.Matrix_A,app.TModel.Assembled.Matrix_B,app.NumberOfEigenValues,app.Shift);SolveTime=toc;
        app.EigenVectors=V;app.EigenValues=diag(D);NewMessage(app,"Solution Time : " + num2str(SolveTime));pause(0.1);
        app.TModel.Solution.EigenValues=app.EigenValues;app.TModel.Solution.EigenVectors=app.EigenVectors;
    end
    app.TModel.Cond=6;
end

            
         

