function [] = ExportAssembledMatrices_EigenModeModule(app)
    filter = {'*.mat'};[file,location] = uiputfile(filter);
    filename=string(location)+string(file);filename=char(filename);
    if(app.TModel.Assembled.Bian),ExportBianisotropicMatrices_EigenMode(filename,app.TModel.Edges,app.TModel.Facets,app.TModel.Assembled.TE,app.TModel.Assembled.TB,app.TModel.Assembled.AM ...
                                                                        ,app.TModel.Assembled.FM,app.TModel.Assembled.FW,app.TModel.Assembled.AW,app.TModel.Assembled.TS,app.TModel.Assembled.TG ...
                                                                        ,app.TModel.Assembled.TBC,app.TModel.Assembled.P,app.TModel.Assembled.TA,app.TModel.Assembled.TC,app.TModel.Assembled.K, ...
                                                                        app.TModel.Assembled.Matrix_A,app.TModel.Assembled.Matrix_B);
    else,ExportMatrices_EigenMode(filename,app.TModel.Edges,app.TModel.Facets,app.TModel.Assembled.TE,app.TModel.Assembled.TB,app.TModel.Assembled.AM,app.TModel.Assembled.FM ...
                                  ,app.TModel.Assembled.FW,app.TModel.Assembled.AW,app.TModel.Assembled.TS,app.TModel.Assembled.TG,app.TModel.Assembled.TBC,app.TModel.Assembled.Matrix_A,app.TModel.Assembled.Matrix_B);
    end
    NewMessage(app,"Model Assembled Matrices Saved in : " + filename);
end

function [] = ExportMatrices_EigenMode(filename,Edges,Facets,TE,TB,A,F,FW,AW,TS,TG,TBC,Matrix_A,Matrix_B)
    save(filename,'Edges','Facets','TE','TB','A','F','FW','AW','TS','TG','TBC','Matrix_A','Matrix_B');
end

function [] = ExportBianisotropicMatrices_EigenMode(filename,Edges,Facets,TE,TB,A,F,FW,AW,TS,TG,TBC,P,TA,TC,K,Matrix_A,Matrix_B)
    save(filename,'Edges','Facets','TE','TB','A','F','FW','AW','TS','TG','TBC','P','TA','TC','K','Matrix_A','Matrix_B');
end