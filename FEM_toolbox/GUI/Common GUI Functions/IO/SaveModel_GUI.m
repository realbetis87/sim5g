function [] = SaveModel_GUI(app)
    app.EigenModeUIFigure.WindowState="minimized";
    filter = {'*.mat'};[file,location] = uiputfile(filter);
    filename=string(location)+string(file);filename=char(filename);
    SaveToolboxModel(filename,app.TModel);
    NewMessage(app,"Model Saved in : " + filename);
end

function [] = SaveToolboxModel(filename,TModel)
    save(filename,'TModel');
end