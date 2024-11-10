%--------------------------------------------------------------------------
%{
    Import ToolboxModel object named TModel

    res = 0 : successfull ToolboxModel Import
    res = 1 : TModel is not a ToolboxModel object
    res = 2 : TModel variable does not exist in filename.mat file
%}
%--------------------------------------------------------------------------
function [TModel,res] = ImportToolboxModel(filename)
    try 
        load(filename,"TModel");
        if(isa(TModel,"ToolboxModel")),res=0;
        else,res=1;
        end
    catch
        res=2;TModel=[];
    end
end

