function [] = BianisotropicCheck_ExcitationModule(app)
            value = app.BianisotropicCheckBox.Value;
            if(value),app.IsotropicCheckBox.Value=false;app.TModel.Domains(app.CurrentDomainIndex).Medium=app.TModel.Domains(app.CurrentDomainIndex).Medium.Bianisotropic();
                app.KsiButton.Enable=true;app.KsiButton.Visible=true;
                app.ZiButton.Enable=true;app.ZiButton.Visible=true;
            else
                app.KsiButton.Enable=false;app.KsiButton.Visible=false;
                app.ZiButton.Enable=false;app.ZiButton.Visible=false;
                if(app.IsotropicCheckBox.Value==false),app.TModel.Domains(app.CurrentDomainIndex).Medium=app.TModel.Domains(app.CurrentDomainIndex).Medium.Anisotropic();
                else,app.TModel.Domains(app.CurrentDomainIndex).Medium=app.TModel.Domains(app.CurrentDomainIndex).Medium.Isotropic();
                end
            end
end

