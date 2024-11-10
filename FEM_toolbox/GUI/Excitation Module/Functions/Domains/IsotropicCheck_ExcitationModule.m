function [] = IsotropicCheck_ExcitationModule(app)
            value = app.IsotropicCheckBox.Value;
            if(value)
                app.BianisotropicCheckBox.Value=false;
                app.TModel.Domains(app.CurrentDomainIndex).Medium=app.TModel.Domains(app.CurrentDomainIndex).Medium.Isotropic();
                app.KsiButton.Enable=false;app.KsiButton.Visible=false;
                app.ZiButton.Enable=false;app.ZiButton.Visible=false;
            else
                if(app.BianisotropicCheckBox.Value==false),app.TModel.Domains(app.CurrentDomainIndex).Medium=app.TModel.Domains(app.CurrentDomainIndex).Medium.Anisotropic();
                else
                end
            end
end

