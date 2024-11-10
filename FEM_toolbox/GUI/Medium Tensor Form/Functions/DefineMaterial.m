function [] = DefineMaterial(varargin)
    if nargin==1,app=varargin{1};
        if(app.medium.Type=="Iso")
            switch app.type
                case 1,app.medium.Epsilon=ReturnFieldNumber(app.YYField);
                case 2,app.medium.Mu=ReturnFieldNumber(app.YYField);
            end
        else
            switch app.type
                case 1
                    app.medium.Epsilon(1,1)=ReturnFieldNumber(app.XXField);
                    app.medium.Epsilon(1,2)=ReturnFieldNumber(app.XYField);
                    app.medium.Epsilon(1,3)=ReturnFieldNumber(app.XZField);
                    app.medium.Epsilon(2,1)=ReturnFieldNumber(app.YXField);
                    app.medium.Epsilon(2,2)=ReturnFieldNumber(app.YYField);
                    app.medium.Epsilon(2,3)=ReturnFieldNumber(app.YZField);
                    app.medium.Epsilon(3,1)=ReturnFieldNumber(app.ZXField);
                    app.medium.Epsilon(3,2)=ReturnFieldNumber(app.ZYField);
                    app.medium.Epsilon(3,3)=ReturnFieldNumber(app.ZZField);
                case 2
                    app.medium.Mu(1,1)=ReturnFieldNumber(app.XXField);
                    app.medium.Mu(1,2)=ReturnFieldNumber(app.XYField);
                    app.medium.Mu(1,3)=ReturnFieldNumber(app.XZField);
                    app.medium.Mu(2,1)=ReturnFieldNumber(app.YXField);
                    app.medium.Mu(2,2)=ReturnFieldNumber(app.YYField);
                    app.medium.Mu(2,3)=ReturnFieldNumber(app.YZField);
                    app.medium.Mu(3,1)=ReturnFieldNumber(app.ZXField);
                    app.medium.Mu(3,2)=ReturnFieldNumber(app.ZYField);
                    app.medium.Mu(3,3)=ReturnFieldNumber(app.ZZField);
                case 3
                    app.medium.Zita(1,1)=ReturnFieldNumber(app.XXField);
                    app.medium.Zita(1,2)=ReturnFieldNumber(app.XYField);
                    app.medium.Zita(1,3)=ReturnFieldNumber(app.XZField);
                    app.medium.Zita(2,1)=ReturnFieldNumber(app.YXField);
                    app.medium.Zita(2,2)=ReturnFieldNumber(app.YYField);
                    app.medium.Zita(2,3)=ReturnFieldNumber(app.YZField);
                    app.medium.Zita(3,1)=ReturnFieldNumber(app.ZXField);
                    app.medium.Zita(3,2)=ReturnFieldNumber(app.ZYField);
                    app.medium.Zita(3,3)=ReturnFieldNumber(app.ZZField);
                case 4
                    app.medium.Ksi(1,1)=ReturnFieldNumber(app.XXField);
                    app.medium.Ksi(1,2)=ReturnFieldNumber(app.XYField);
                    app.medium.Ksi(1,3)=ReturnFieldNumber(app.XZField);
                    app.medium.Ksi(2,1)=ReturnFieldNumber(app.YXField);
                    app.medium.Ksi(2,2)=ReturnFieldNumber(app.YYField);
                    app.medium.Ksi(2,3)=ReturnFieldNumber(app.YZField);
                    app.medium.Ksi(3,1)=ReturnFieldNumber(app.ZXField);
                    app.medium.Ksi(3,2)=ReturnFieldNumber(app.ZYField);
                    app.medium.Ksi(3,3)=ReturnFieldNumber(app.ZZField);
            end
        end
    elseif nargin==2,app=varargin{1};index=varargin{2};
        if(app.medium.Type=="Iso")
            switch app.type
                case 1,app.medium.Epsilon(index)=ReturnFieldNumber(app.YYField);
                case 2,app.medium.Mu(index)=ReturnFieldNumber(app.YYField);
            end
        else
            switch app.type
                case 1,er=zeros(3,3);
                       er(1,1)=ReturnFieldNumber(app.XXField);
                       er(1,2)=ReturnFieldNumber(app.XYField);
                       er(1,3)=ReturnFieldNumber(app.XZField);
                       er(2,1)=ReturnFieldNumber(app.YXField);
                       er(2,2)=ReturnFieldNumber(app.YYField);
                       er(2,3)=ReturnFieldNumber(app.YZField);
                       er(3,1)=ReturnFieldNumber(app.ZXField);
                       er(3,2)=ReturnFieldNumber(app.ZYField);
                       er(3,3)=ReturnFieldNumber(app.ZZField);
                       app.medium.Epsilon{index}=er;
                case 2,mr=zeros(3,3);
                       mr(1,1)=ReturnFieldNumber(app.XXField);
                       mr(1,2)=ReturnFieldNumber(app.XYField);
                       mr(1,3)=ReturnFieldNumber(app.XZField);
                       mr(2,1)=ReturnFieldNumber(app.YXField);
                       mr(2,2)=ReturnFieldNumber(app.YYField);
                       mr(2,3)=ReturnFieldNumber(app.YZField);
                       mr(3,1)=ReturnFieldNumber(app.ZXField);
                       mr(3,2)=ReturnFieldNumber(app.ZYField);
                       mr(3,3)=ReturnFieldNumber(app.ZZField);
                       app.medium.Mu{index}=mr;
                case 3,zita=zeros(3,3);
                       zita(1,1)=ReturnFieldNumber(app.XXField);
                       zita(1,2)=ReturnFieldNumber(app.XYField);
                       zita(1,3)=ReturnFieldNumber(app.XZField);
                       zita(2,1)=ReturnFieldNumber(app.YXField);
                       zita(2,2)=ReturnFieldNumber(app.YYField);
                       zita(2,3)=ReturnFieldNumber(app.YZField);
                       zita(3,1)=ReturnFieldNumber(app.ZXField);
                       zita(3,2)=ReturnFieldNumber(app.ZYField);
                       zita(3,3)=ReturnFieldNumber(app.ZZField);
                       app.medium.Zita{index}=zita;
                case 4,ksi=zeros(3,3);
                       ksi(1,1)=ReturnFieldNumber(app.XXField);
                       ksi(1,2)=ReturnFieldNumber(app.XYField);
                       ksi(1,3)=ReturnFieldNumber(app.XZField);
                       ksi(2,1)=ReturnFieldNumber(app.YXField);
                       ksi(2,2)=ReturnFieldNumber(app.YYField);
                       ksi(2,3)=ReturnFieldNumber(app.YZField);
                       ksi(3,1)=ReturnFieldNumber(app.ZXField);
                       ksi(3,2)=ReturnFieldNumber(app.ZYField);
                       ksi(3,3)=ReturnFieldNumber(app.ZZField);
                       app.medium.Ksi{index}=ksi;
            end
        end

    end
end

function [res] = ReturnFieldNumber(field),res=String2Complex(field.Value);end%res=str2double(res);end

