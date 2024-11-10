function [] = SaveFields_ExcitationMode(app,Mode)
        app.ExcitationUIFigure.WindowState="minimized";
        filter = {'*.mat'};[file,location] = uiputfile(filter);
        filename=string(location)+string(file);filename=char(filename);
        if(Mode==0)
            if(app.TModel.Frequency.NF==1)
                X=app.TModel.Solution.SolutionVector;
                Nv=get(app.HNField,"Value");Nh=get(app.VNField,"Value");Pos=get(app.PositionField,"Value");
                AxisIndex=app.PlotAxisSelection.ValueIndex;
                switch AxisIndex
                    case 1,position=[Pos;nan;nan;];
                    case 2,position=[nan;Pos;nan;];
                    case 3,position=[nan;nan;Pos;];
                end
                [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz] = ReturnFields_Excitation(app.TModel,X,position,[Nv,Nh]);
                SaveFields(filename,Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz);
                NewMessage(app,"Fields Saved in : " + filename);
            else
                freq_index=app.ResFreq.ValueIndex;X=app.TModel.Solution.SolutionVector(:,freq_index);
                Nv=get(app.HNField,"Value");Nh=get(app.VNField,"Value");Pos=get(app.PositionField,"Value");
                AxisIndex=app.PlotAxisSelection.ValueIndex;
                switch AxisIndex
                       case 1,position=[Pos;nan;nan;];
                       case 2,position=[nan;Pos;nan;];
                       case 3,position=[nan;nan;Pos;];
                end
                [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz] = ReturnFields_Excitation(app.TModel,X,position,[Nv,Nh]);
                SaveFields(filename,Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz);
                NewMessage(app,"Fields Saved in : " + filename);
            end
        else
               SaveAll(filename,app.TModel.Solution.SolutionVector,app.TModel);
        end
        app.ExcitationUIFigure.WindowState="normal";
end


function [] = SaveFields(filename,Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz)
    save(filename,'Ex','Ey','Ez','Bx','By','Bz','xx','yy','zz');
end
function [] = SaveAll(filename,SolutionVector,toolboxModel)
save(filename,'SolutionVector','toolboxModel');
end