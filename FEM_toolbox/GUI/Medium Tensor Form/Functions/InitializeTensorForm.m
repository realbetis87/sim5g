function [] = InitializeTensorForm(app),medium=app.medium;app.NameField.Value=app.medium.Tag;
    %----------------------------------------------------------------------
    switch app.type
        case 1,app.TitleLabel.Text="Dielectric Permittivity ε_{r}";app.ExplainLabel.Text="";
        case 2,app.TitleLabel.Text="Magnetic Permability μ_{r}";app.ExplainLabel.Text="";
        case 3,app.TitleLabel.Text="Magneto electric coupling ζ";app.ExplainLabel.Text="\textbf{Β}=μ_{0}μ_{r}\textbf{Η}+ζ\textbf{Ε}";
        case 4,app.TitleLabel.Text="Magneto electric coupling ξ";app.ExplainLabel.Text="\textbf{D}=ε_{0}ε_{r}\textbf{E}+ξ\textbf{H}";
     end
    %----------------------------------------------------------------------
    if(~isempty(app.masterApp.TModel.Frequency)),app.MultiF=app.masterApp.TModel.Frequency;end
    %----------------------------------------------------------------------
    if(app.medium.IsDispersive)
        switch medium.FRange.Unit
            case "Hz",mul=1;
            case "KHz",mul=1e-3;
            case "MHz",mul=1e-6;
            case "GHz",mul=1e-9;
            case "THz",mul=1e-12;
        end
        for ii=1:medium.FRange.NF,app.DropDown.Items{ii}=num2str(medium.FRange.Frequency(ii)*mul);end,app.UnitLabel.Text=medium.FRange.Unit;app.DropDown.ValueIndex=1;
            if(medium.Type=="Iso")
                app.RXLabel.Visible=false;app.RYLabel.Visible=false;app.RZLabel.Visible=false;
                app.CXLabel.Visible=false;app.CYLabel.Visible=false;app.CZLabel.Visible=false;
                app.XXField.Enable=false;app.XYField.Enable=false;app.XZField.Enable=false;
                app.ZXField.Enable=false;app.ZYField.Enable=false;app.ZZField.Enable=false;
                app.YXField.Enable=false;app.YZField.Enable=false;
                switch app.type
                    case 1,set(app.YYField,"Value",num2str(medium.Epsilon(1)));
                    case 2,set(app.YYField,"Value",num2str(medium.Mu(1)));
                end
            else
                switch app.type
                    case 1,epsilon=medium.Epsilon{1};
                        set(app.XXField,"Value",num2str(epsilon(1,1)));set(app.XYField,"Value",num2str(epsilon(1,2)));set(app.XZField,"Value",num2str(epsilon(1,3)));
                        set(app.YXField,"Value",num2str(epsilon(2,1)));set(app.YYField,"Value",num2str(epsilon(2,2)));set(app.YZField,"Value",num2str(epsilon(2,3)));
                        set(app.ZXField,"Value",num2str(epsilon(3,1)));set(app.ZYField,"Value",num2str(epsilon(3,2)));set(app.ZZField,"Value",num2str(epsilon(3,3)));
                    case 2,mu=medium.Mu{1};
                        set(app.XXField,"Value",num2str(mu(1,1)));set(app.XYField,"Value",num2str(mu(1,2)));set(app.XZField,"Value",num2str(mu(1,3)));
                        set(app.YXField,"Value",num2str(mu(2,1)));set(app.YYField,"Value",num2str(mu(2,2)));set(app.YZField,"Value",num2str(mu(2,3)));
                        set(app.ZXField,"Value",num2str(mu(3,1)));set(app.ZYField,"Value",num2str(mu(3,2)));set(app.ZZField,"Value",num2str(mu(3,3)));
                    case 3,zita=medium.Zita{1};
                        set(app.XXField,"Value",num2str(zita(1,1)));set(app.XYField,"Value",num2str(zita(1,2)));set(app.XZField,"Value",num2str(zita(1,3)));
                        set(app.YXField,"Value",num2str(zita(2,1)));set(app.YYField,"Value",num2str(zita(2,2)));set(app.YZField,"Value",num2str(zita(2,3)));
                        set(app.ZXField,"Value",num2str(zita(3,1)));set(app.ZYField,"Value",num2str(zita(3,2)));set(app.ZZField,"Value",num2str(zita(3,3)));
                    case 4,ksi=medium.Ksi{1};
                        set(app.XXField,"Value",num2str(ksi(1,1)));set(app.XYField,"Value",num2str(ksi(1,2)));set(app.XZField,"Value",num2str(ksi(1,3)));
                        set(app.YXField,"Value",num2str(ksi(2,1)));set(app.YYField,"Value",num2str(ksi(2,2)));set(app.YZField,"Value",num2str(ksi(2,3)));
                        set(app.ZXField,"Value",num2str(ksi(3,1)));set(app.ZYField,"Value",num2str(ksi(3,2)));set(app.ZZField,"Value",num2str(ksi(3,3)));
                end
            end
    else
        app.MFPanel.Visible=false;app.MFPanel.Enable=false;
        if(medium.Type=="Iso")
            app.RXLabel.Visible=false;app.RYLabel.Visible=false;app.RZLabel.Visible=false;
            app.CXLabel.Visible=false;app.CYLabel.Visible=false;app.CZLabel.Visible=false;
            app.XXField.Enable=false;app.XYField.Enable=false;app.XZField.Enable=false;
            app.ZXField.Enable=false;app.ZYField.Enable=false;app.ZZField.Enable=false;
            app.YXField.Enable=false;app.YZField.Enable=false;
            %app.TensorPanel.RowHeight={"0x","0x","1x","0x"};app.TensorPanel.ColumnWidth={"0x","0x","1x","0x"};
            switch app.type
                case 1,set(app.YYField,"Value",num2str(medium.Epsilon));
                case 2,set(app.YYField,"Value",num2str(medium.Mu));
            end
        else
            switch app.type
                    case 1,epsilon=medium.Epsilon;
                        set(app.XXField,"Value",num2str(epsilon(1,1)));set(app.XYField,"Value",num2str(epsilon(1,2)));set(app.XZField,"Value",num2str(epsilon(1,3)));
                        set(app.YXField,"Value",num2str(epsilon(2,1)));set(app.YYField,"Value",num2str(epsilon(2,2)));set(app.YZField,"Value",num2str(epsilon(2,3)));
                        set(app.ZXField,"Value",num2str(epsilon(3,1)));set(app.ZYField,"Value",num2str(epsilon(3,2)));set(app.ZZField,"Value",num2str(epsilon(3,3)));
                    case 2,mu=medium.Mu;
                        set(app.XXField,"Value",num2str(mu(1,1)));set(app.XYField,"Value",num2str(mu(1,2)));set(app.XZField,"Value",num2str(mu(1,3)));
                        set(app.YXField,"Value",num2str(mu(2,1)));set(app.YYField,"Value",num2str(mu(2,2)));set(app.YZField,"Value",num2str(mu(2,3)));
                        set(app.ZXField,"Value",num2str(mu(3,1)));set(app.ZYField,"Value",num2str(mu(3,2)));set(app.ZZField,"Value",num2str(mu(3,3)));
                    case 3,zita=medium.Zita;
                        set(app.XXField,"Value",num2str(zita(1,1)));set(app.XYField,"Value",num2str(zita(1,2)));set(app.XZField,"Value",num2str(zita(1,3)));
                        set(app.YXField,"Value",num2str(zita(2,1)));set(app.YYField,"Value",num2str(zita(2,2)));set(app.YZField,"Value",num2str(zita(2,3)));
                        set(app.ZXField,"Value",num2str(zita(3,1)));set(app.ZYField,"Value",num2str(zita(3,2)));set(app.ZZField,"Value",num2str(zita(3,3)));
                    case 4,ksi=medium.Ksi;
                        set(app.XXField,"Value",num2str(ksi(1,1)));set(app.XYField,"Value",num2str(ksi(1,2)));set(app.XZField,"Value",num2str(ksi(1,3)));
                        set(app.YXField,"Value",num2str(ksi(2,1)));set(app.YYField,"Value",num2str(ksi(2,2)));set(app.YZField,"Value",num2str(ksi(2,3)));
                        set(app.ZXField,"Value",num2str(ksi(3,1)));set(app.ZYField,"Value",num2str(ksi(3,2)));set(app.ZZField,"Value",num2str(ksi(3,3)));
           end
        end
    end
end

