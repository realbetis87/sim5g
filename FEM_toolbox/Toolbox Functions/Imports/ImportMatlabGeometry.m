function [TModel,res] = ImportMatlabGeometry(filename)
    try
        load(filename,'model');
    catch
    end
end

