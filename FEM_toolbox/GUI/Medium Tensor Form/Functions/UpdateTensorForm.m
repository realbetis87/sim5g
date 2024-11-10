function [] = UpdateTensorForm(app,index),medium=app.medium;
        if(medium.Type=="Iso")
            switch app.type
                case 1,set(app.YYField,"Value",num2str(medium.Epsilon(index)));
                case 2,set(app.YYField,"Value",num2str(medium.Mu(index)));
            end
        else
            switch app.type
                case 1,epsilon=medium.Epsilon{index};
                        set(app.XXField,"Value",num2str(epsilon(1,1)));set(app.XYField,"Value",num2str(epsilon(1,2)));set(app.XZField,"Value",num2str(epsilon(1,3)));
                        set(app.YXField,"Value",num2str(epsilon(2,1)));set(app.YYField,"Value",num2str(epsilon(2,2)));set(app.YZField,"Value",num2str(epsilon(2,3)));
                        set(app.ZXField,"Value",num2str(epsilon(3,1)));set(app.ZYField,"Value",num2str(epsilon(3,2)));set(app.ZZField,"Value",num2str(epsilon(3,3)));
                case 2,mu=medium.Mu{index};
                        set(app.XXField,"Value",num2str(mu(1,1)));set(app.XYField,"Value",num2str(mu(1,2)));set(app.XZField,"Value",num2str(mu(1,3)));
                        set(app.YXField,"Value",num2str(mu(2,1)));set(app.YYField,"Value",num2str(mu(2,2)));set(app.YZField,"Value",num2str(mu(2,3)));
                        set(app.ZXField,"Value",num2str(mu(3,1)));set(app.ZYField,"Value",num2str(mu(3,2)));set(app.ZZField,"Value",num2str(mu(3,3)));
                case 3,zita=medium.Zita{index};
                        set(app.XXField,"Value",num2str(zita(1,1)));set(app.XYField,"Value",num2str(zita(1,2)));set(app.XZField,"Value",num2str(zita(1,3)));
                        set(app.YXField,"Value",num2str(zita(2,1)));set(app.YYField,"Value",num2str(zita(2,2)));set(app.YZField,"Value",num2str(zita(2,3)));
                        set(app.ZXField,"Value",num2str(zita(3,1)));set(app.ZYField,"Value",num2str(zita(3,2)));set(app.ZZField,"Value",num2str(zita(3,3)));
                case 4,ksi=medium.Ksi{index};
                        set(app.XXField,"Value",num2str(ksi(1,1)));set(app.XYField,"Value",num2str(ksi(1,2)));set(app.XZField,"Value",num2str(ksi(1,3)));
                        set(app.YXField,"Value",num2str(ksi(2,1)));set(app.YYField,"Value",num2str(ksi(2,2)));set(app.YZField,"Value",num2str(ksi(2,3)));
                        set(app.ZXField,"Value",num2str(ksi(3,1)));set(app.ZYField,"Value",num2str(ksi(3,2)));set(app.ZZField,"Value",num2str(ksi(3,3)));
            end
        end
end

