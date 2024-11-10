function [] = ExportAssembledMatrices_ExcitationModule(app)
    filter = {'*.mat'};[file,location] = uiputfile(filter);
    filename=string(location)+string(file);filename=char(filename);
    if(app.TModel.Assembled.Bian)
        if(app.TModel.Assembled.IsDir)
                                        ExportBianisotropicMatrices_Excitation(filename,app.TModel.Edges,app.TModel.Facets,app.TModel.Assembled.TE,app.TModel.Assembled.TB,app.TModel.Assembled.AM ...
                                                                        ,app.TModel.Assembled.FM,app.TModel.Assembled.TP,app.TModel.Assembled.TPV,app.TModel.Assembled.TS,app.TModel.Assembled.TG ...
                                                                        ,app.TModel.Assembled.TBC,app.TModel.Assembled.P,app.TModel.Assembled.TA,app.TModel.Assembled.TC ...
                                                                        ,app.TModel.Assembled.TEB,app.TModel.Assembled.TBB,app.TModel.Assembled.FB,app.TModel.Assembled.AB...
                                                                        ,app.TModel.Assembled.PB,app.TModel.Assembled.TAB,app.TModel.Assembled.TCB,...
                                                                        app.TModel.Assembled.Matrix_A,app.TModel.Assembled.Matrix_B);
        else
                                        ExportBianisotropicMatrices_DirMode(filename,app.TModel.Edges,app.TModel.Facets,app.TModel.Assembled.TE,app.TModel.Assembled.TB,app.TModel.Assembled.AM ...
                                                                        ,app.TModel.Assembled.FM,app.TModel.Assembled.TP,app.TModel.Assembled.TPV,app.TModel.Assembled.TS,app.TModel.Assembled.TG ...
                                                                        ,app.TModel.Assembled.TBC,app.TModel.Assembled.P,app.TModel.Assembled.TA,app.TModel.Assembled.TC, ...
                                                                        app.TModel.Assembled.Matrix_A,app.TModel.Assembled.Matrix_B);
        end
    else
        if(app.TModel.Assembled.IsDir)
                                        ExportMatrices_Excitation(filename,app.TModel.Edges,app.TModel.Facets,app.TModel.Assembled.TE,app.TModel.Assembled.TB,app.TModel.Assembled.AM,app.TModel.Assembled.FM, ...
                                        app.TModel.Assembled.TP,app.TModel.Assembled.TPV,app.TModel.Assembled.TS,app.TModel.Assembled.TG,app.TModel.Assembled.TBC,app.TModel.Assembled.Matrix_A,app.TModel.Assembled.Matrix_B);
        else
        end
                                        ExportMatrices_DirExcitation(filename,app.TModel.Edges,app.TModel.Facets,app.TModel.Assembled.TE,app.TModel.Assembled.TB,app.TModel.Assembled.AM,app.TModel.Assembled.FM ...
                                            ,app.TModel.Assembled.TP,app.TModel.Assembled.TPV,app.TModel.Assembled.TS,app.TModel.Assembled.TG,app.TModel.Assembled.TBC ...
                                            ,app.TModel.Assembled.TEB,app.TModel.Assembled.TBB,app.TModel.Assembled.AB,app.TModel.Assembled.FB,app.TModel.Assembled.Matrix_A,app.TModel.Assembled.Matrix_B);
    end
    NewMessage(app,"Model Assembled Matrices Saved in : " + filename);
end

function [] = ExportMatrices_Excitation(filename,Edges,Facets,TE,TB,A,F,TP,TPV,TS,TG,TBC,Matrix_A,Matrix_B)
    save(filename,'Edges','Facets','TE','TB','A','F','TP','TPV','TS','TG','TBC','Matrix_A','Matrix_B');
end
function [] = ExportMatrices_DirExcitation(filename,Edges,Facets,TE,TB,A,F,TP,TPV,TS,TG,TBC,TEB,TBB,AB,FB,Matrix_A,Matrix_B)
    save(filename,'Edges','Facets','TE','TB','A','F','TP','TPV','TS','TG','TBC','TEB','TBB','AB','FB','Matrix_A','Matrix_B');
end

function [] = ExportBianisotropicMatrices_Excitation(filename,Edges,Facets,TE,TB,A,F,TP,TPV,TS,TG,TBC,P,TA,TC,Matrix_A,Matrix_B)
    save(filename,'Edges','Facets','TE','TB','A','F','TP','TPV','TS','TG','TBC','P','TA','TC','Matrix_A','Matrix_B');
end

function [] = ExportBianisotropicMatrices_DirMode(filename,Edges,Facets,TE,TB,A,F,TP,TPV,TS,TG,TBC,P,TA,TC,K,TEB,TBB,AB,FB,PB,TAB,TCB,Matrix_A,Matrix_B)
    save(filename,'Edges','Facets','TE','TB','A','F','TP','TPV','TS','TG','TBC','P','TA','TC','K','TEB','TBB','AB','FB','AB','PB','TAB',"TCB",'Matrix_A','Matrix_B');
end