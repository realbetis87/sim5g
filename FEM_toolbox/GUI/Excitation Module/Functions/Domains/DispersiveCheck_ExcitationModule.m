function [] = DispersiveCheck_ExcitationModule(app)
     value = app.DispersiveCheckBox.Value;
     if(value),app.TModel.Domains(app.CurrentDomainIndex).Medium=app.TModel.Domains(app.CurrentDomainIndex).Medium.Dispersive(app.TModel.Frequency);
     else,app.TModel.Domains(app.CurrentDomainIndex).Medium=app.TModel.Domains(app.CurrentDomainIndex).Medium.NoDispersive();
     end  
end

