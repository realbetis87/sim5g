function [] = GetBoundaryTensors(app)
            switch app.mode
                case 0, app.boundary.Param=ReturnFieldNumber(app.TensorField);
                case 1,app.boundary.Param(app.currentFrequencyIndex)=ReturnFieldNumber(app.TensorField);
                case 2,tensor=zeros(3,3);tensor(1,1)=ReturnFieldNumber(app.XXField);tensor(1,2)=ReturnFieldNumber(app.XYField);tensor(1,3)=ReturnFieldNumber(app.XZField);
                                         tensor(2,1)=ReturnFieldNumber(app.YXField);tensor(2,2)=ReturnFieldNumber(app.YYField);tensor(2,3)=ReturnFieldNumber(app.YZField);
                                         tensor(3,1)=ReturnFieldNumber(app.ZXField);tensor(3,2)=ReturnFieldNumber(app.ZYField);tensor(3,3)=ReturnFieldNumber(app.ZZField);
                       app.boundary.Param=tensor;
                case 3,tensor=zeros(3,3);tensor(1,1)=ReturnFieldNumber(app.XXField);tensor(1,2)=ReturnFieldNumber(app.XYField);tensor(1,3)=ReturnFieldNumber(app.XZField);
                                         tensor(2,1)=ReturnFieldNumber(app.YXField);tensor(2,2)=ReturnFieldNumber(app.YYField);tensor(2,3)=ReturnFieldNumber(app.YZField);
                                         tensor(3,1)=ReturnFieldNumber(app.ZXField);tensor(3,2)=ReturnFieldNumber(app.ZYField);tensor(3,3)=ReturnFieldNumber(app.ZZField);
                       app.boundary.Param{app.currentFrequencyIndex}=tensor;
                case 4,app.boundary.Param=str2double(app.SlaveSelection.Value);
            end
end
function [res] = ReturnFieldNumber(field),res=String2Complex(field.Value);end