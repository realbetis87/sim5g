function [] = Solve_2DModalSolver(app)
    if (app.TModel.Frequency.NF==1)
            app.NumberOfEigenValues=get(app.NEigenValueField,"Value");NE=app.NumberOfEigenValues;
            app.EigenValueShift=get(app.ShiftField,"Value");shift=app.EigenValueShift;
            M1=app.Excitation_Res.Assembly.Matrix_A;M2=app.Excitation_Res.Assembly.Matrix_B;
            if(app.BetaCheck.Value)
                tic;[app.EigenVectors,app.EigenValues]=eigs(M1,M2,NE,shift);SolTime=toc;
            else,ElectromagneticConstants;k0=2*pi*app.TModel.Frequency.Frequency/c0;
                tic;[app.EigenVectors,app.EigenValues]=eigs(M1,k0*M2,NE,shift);SolTime=toc;
            end
            app.EigenValues=diag(app.EigenValues);
            NewMessage(app,"Solution Time: " + num2str(SolTime));pause(0.1);
    else
        app.EigenValues=cell(app.TModel.Frequency.NF,1);app.EigenVectors=cell(app.TModel.Frequency.NF,1);
        if(app.SolFreq.ValueIndex==1),tic;
            for ii=1:app.TModel.Frequency.NF
                app.NumberOfEigenValues=get(app.NEigenValueField,"Value");NE=app.NumberOfEigenValues;
                app.EigenValueShift=get(app.ShiftField,"Value");shift=app.EigenValueShift;
                M1=app.Excitation_Res.Assembly.Matrix_A{ii};M2=app.Excitation_Res.Assembly.Matrix_B{ii};
                if(app.BetaCheck.Value)
                    [eigenVectors,eigenValues]=eigs(M1,M2,NE,shift);
                else,ElectromagneticConstants;k0=2*pi*app.TModel.Frequency.Frequency(ii)/c0;
                    [eigenVectors,eigenValues]=eigs(M1,k0*M2,NE,shift);
                end
                app.EigenVectors{ii}=eigenVectors;app.EigenValues{ii}=diag(eigenValues);
            end,SolTime=toc;
            NewMessage(app,"Solution Time: " + num2str(SolTime));pause(0.1);
        else
            ii=app.SolFreq.ValueIndex-1;
            app.NumberOfEigenValues=get(app.NEigenValueField,"Value");NE=app.NumberOfEigenValues;
            app.EigenValueShift=get(app.ShiftField,"Value");shift=app.EigenValueShift;
            M1=app.Excitation_Res.Assembly.Matrix_A{ii};M2=app.Excitation_Res.Assembly.Matrix_B{ii};
                if(app.BetaCheck.Value)
                    tic;[eigenVectors,eigenValues]=eigs(M1,M2,NE,shift);SolTime=toc;
                else,ElectromagneticConstants;k0=2*pi*app.TModel.Frequency.Frequency(ii)/c0;
                    tic;[eigenVectors,eigenValues]=eigs(M1,k0*M2,NE,shift);SolTime=toc;
                end
            app.EigenVectors{ii}=eigenVectors;app.EigenValues{ii}=diag(eigenValues);
            NewMessage(app,"Solution Time for Frequency"+ num2str(app.TModel.Frequency.UFrequency)+ app.TModel.Frequency.Unit+ " : " + num2str(SolTime) );pause(0.1);
        end
    end
end

