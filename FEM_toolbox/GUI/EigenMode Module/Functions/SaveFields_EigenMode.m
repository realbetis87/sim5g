function [] = SaveFields_EigenMode(app,Mode)
    if(app.TModel.Frequency.NF==1)
        app.EigenModeUIFigure.WindowState="minimized";
        filter = {'*.mat'};[file,location] = uiputfile(filter);
        filename=string(location)+string(file);filename=char(filename);
        if(Mode==0)
            eig_index=app.ResSelection.ValueIndex;X=app.TModel.Solution.EigenVectors(:,eig_index);eigenValue=app.TModel.Solution.EigenValues(eig_index);
            Nv=get(app.HNField,"Value");Nh=get(app.VNField,"Value");Pos=get(app.PositionField,"Value");
            AxisIndex=app.PlotAxisSelection.ValueIndex;
            switch AxisIndex
                case 1,position=[Pos;nan;nan;];
                case 2,position=[nan;Pos;nan;];
                case 3,position=[nan;nan;Pos;];
            end
            [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz] = ReturnFields_EigenMode(app.TModel,X,position,[Nv;Nh;]);
            SaveFields(filename,eigenValue,Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz)
        else
        SaveAll(filename,app.TModel.Solution.EigenVectors,app.TModel.Solution.EigenValues);
        end
        NewMessage(app,"Fields Saved in : " + filename);
    else
        app.EigenModeUIFigure.WindowState="minimized";
        filter = {'*.mat'};[file,location] = uiputfile(filter);
        filename=string(location)+string(file);filename=char(filename);
        if(Mode==0)
            freq_index=app.ResFreq.ValueIndex;eig_index=app.ResSelection.ValueIndex;
            eigenVectors=app.TModel.Solution.EigenVectors{freq_index};eigenValues=app.TModel.Solution.EigenValues{freq_index};
            X=eigenVectors(:,eig_index);eigenValue=eigenValues(eig_index);
            Nv=get(app.HNField,"Value");Nh=get(app.VNField,"Value");Pos=get(app.PositionField,"Value");
            AxisIndex=app.PlotAxisSelection.ValueIndex;
            switch AxisIndex
                case 1,position=[Pos;nan;nan;];
                case 2,position=[nan;Pos;nan;];
                case 3,position=[nan;nan;Pos;];
            end
            [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz] = ReturnFields_EigenMode(app.TModel,X,position,[Nv;Nh;]);
            SaveFields(filename,eigenValue,Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz)
        else
        SaveAll(filename,app.TModel.Solution.EigenVectors,app.TModel.Solution.EigenValues);
        end
        NewMessage(app,"Fields Saved in : " + filename);
    end
    app.EigenModeUIFigure.WindowState="normal";
end

function [] = SaveFields(filename,eigenValue,Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz)
    save(filename,'eigenValue','Ex','Ey','Ez','Bx','By','Bz','xx','yy','zz');
end
function [] = SaveAll(filename,EigenVectors,EigenValues)
save(filename,'EigenVectors','EigenValues');
end

