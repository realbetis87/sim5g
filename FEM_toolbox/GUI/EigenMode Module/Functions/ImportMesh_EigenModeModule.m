function [] = ImportMesh_EigenModeModule(app)
    if(app.MatlabMeshButton.Value)
         filename=uigetfile('*.mat',"Load Matlab Mesh file");[TModel,res]=ImportMatlabGeometry(filename);
         switch res
             case 0,app.TModel=TModel;app.GeometryLoaded=true;NewMessage(app,"Matlab Mesh Loaded");GeometryLoaded_EigenModeModule(app);GUIActive(app,"Successful Import");app.TModel.Cond=1;
             case 1,NewMessage(app,"Matlab Mesh Import failed.");NewMessage(app,"model variable not a PDEModel Object");GUIActive(app,"Mesh Import Failed");
             case 2,NewMessage(app,"Matlab Mesh Import failed.");NewMessage(app,"model variable not in file");GUIActive(app,"Mesh Import Failed");
         end
    elseif(app.ComsolMeshButton.Value)
         filename=uigetfile('*.mat',"Load Comsol Mesh file");[TModel,res]=ImportComsolGeometry(filename);
         switch res
             case 0,app.TModel=TModel;app.GeometryLoaded=true;NewMessage(app,"Comsol Mesh Loaded"); GeometryLoaded_EigenModeModule(app);GUIActive(app,"Successful Import");app.TModel.Cond=1;
             case 1,NewMessage(app,"Comsol Mesh Import failed.");NewMessage(app,"Check p,t,meshdata variables.");GUIActive(app,"Mesh Import Failed");
             case 2,NewMessage(app,"Comsol Mesh Import failed.");NewMessage(app,"p,t or meshdata variables do not exist in file.");GUIActive(app,"Mesh Import Failed");
         end
    elseif(app.ToolboxModelButton.Value)
         filename=uigetfile('*.mat',"Load Toolbox Model file");[TModel,res]=ImportToolboxModel(filename);
         switch res
             case 0,app.TModel=TModel;app.InitializationsComplete=true;NewMessage(app,"Toolbox Model Loaded");ToolboxModelLoaded_EigenModeModule(app);GUIActive(app,"Successful Import");
             case 1,NewMessage(app,"Toolbox Model Import failed.");NewMessage(app,"TModel variable is not a a ToolboxModel Object");GUIActive(app,"Toolbox Model Import Failed");
             case 2,NewMessage(app,"Toolbox Model Import failed.");NewMessage(app,"TModel variabledoes not exist in file");GUIActive(app,"Toolbox Model Import Failed");
         end
    end
end