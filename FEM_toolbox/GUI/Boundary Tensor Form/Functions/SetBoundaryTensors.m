function [] = SetBoundaryTensors(app)
            switch app.mode
                case 0,set(app.TensorField,"Value",num2str(app.boundary.Param));
                case 1,set(app.TensorField,"Value", num2str(app.boundary.Param(app.currentFrequencyIndex)));
                case 2,set(app.XXField,"Value",num2str(app.boundary.Param(1,1)));set(app.XYField,"Value",num2str(app.boundary.Param(1,2)));set(app.XZField,"Value",num2str(app.boundary.Param(1,3)));
                       set(app.YXField,"Value",num2str(app.boundary.Param(2,1)));set(app.YYField,"Value",num2str(app.boundary.Param(2,2)));set(app.YZField,"Value",num2str(app.boundary.Param(2,3)));
                       set(app.ZXField,"Value",num2str(app.boundary.Param(3,1)));set(app.ZYField,"Value",num2str(app.boundary.Param(3,2)));set(app.ZZField,"Value",num2str(app.boundary.Param(3,3)));
                case 3,tensor=app.boundary.Param{app.currentFrequencyIndex};set(app.XXField,"Value",num2str(tensor(1,1)));set(app.XYField,"Value",num2str(tensor(1,2)));set(app.XZField,"Value",num2str(tensor(1,3)));
                                                                            set(app.YXField,"Value",num2str(tensor(2,1)));set(app.YYField,"Value",num2str(tensor(2,2)));set(app.YZField,"Value",num2str(tensor(2,3)));
                                                                            set(app.ZXField,"Value",num2str(tensor(3,1)));set(app.ZYField,"Value",num2str(tensor(3,2)));set(app.ZZField,"Value",num2str(tensor(3,3)));
                       
            end
end

