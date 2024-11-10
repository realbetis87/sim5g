function [] = UpdateCurrentDomain(app),PlotDomain_GUI(app,app.CurrentDomainIndex);
    medium=app.TModel.Domains(app.CurrentDomainIndex).Medium;app.MediumLabel.Text=medium.Tag;
    if(medium.Type=="Iso"),app.IsotropicCheckBox.Value=true;app.BianisotropicCheckBox.Value=false;app.KsiButton.Enable=false;app.KsiButton.Visible=false;app.ZiButton.Enable=false;app.ZiButton.Visible=false;
    elseif(medium.Type=="Anis"),app.IsotropicCheckBox.Value=false;app.BianisotropicCheckBox.Value=false;app.KsiButton.Enable=false;app.KsiButton.Visible=false;app.ZiButton.Enable=false;app.ZiButton.Visible=false;
    elseif(medium.Type=="Bian"),app.IsotropicCheckBox.Value=false;app.BianisotropicCheckBox.Value=true;app.KsiButton.Enable=true;app.KsiButton.Visible=true;app.ZiButton.Enable=true;app.ZiButton.Visible=true;
    end
    if(isempty(medium.IsDispersive)),app.DispersiveCheckBox.Value=false;
    elseif(medium.IsDispersive==false),app.DispersiveCheckBox.Value=false;
    elseif(medium.IsDispersive),app.DispersiveCheckBox.Value=true;
    end
end
