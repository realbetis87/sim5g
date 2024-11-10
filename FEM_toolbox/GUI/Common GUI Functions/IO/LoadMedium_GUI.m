function [] = LoadMedium_GUI(app)
     filename=uigetfile('*.mat',"Load Medium file");
     [res] = LoadMedium(app,filename);
     UpdateCurrentDomain(app);
     NewMessage(app,res);
end


function [res] = LoadMedium(app,filename)
    try
        load(filename,'medium');
        if(isa(medium,"Medium"))
            app.TModel.Domains(app.CurrentDomainIndex).Medium=medium;
            res = "Domain " + num2str(app.CurrentDomainIndex) + " medium loaded";
        else
          res = "Medium load failed, medium variable is not a toolbox object"; 
        end
    catch
         try 
             load(filename,'epsilon','mu','ksi','zita','tag');
             medium=app.TModel.Domains(app.CurrentDomainIndex).Medium;medium.Tag=tag;
             if(ismatrix(epsilon) && ismatrix(mu) && ismatrix(ksi) && ismatrix(zita))
                if(size(epsilon,1)==3 && size(epsilon,2)==3 && size(mu,1)==3 && size(mu,2)==3 && size(ksi,1)==3 && size(ksi,2)==3 && size(zita,1)==3 && size(zita,2)==3)
                   medium=medium.NoDispersive();medium=medium.Bianisotropic();
                   medium.Epsilon=epsilon;medium.Mu=mu;medium.Ksi=ksi;medium.Zita=zita;
                   app.TModel.Domains(app.CurrentDomainIndex).Medium=medium;
                   res = "Non Dispersive Bianisotropic Medium  for Domain "+ num2str(app.CurrentDomainIndex) +" Sucessfully loaded";
                else
                        res=" Load Unsuccessful Tensor Matrices not 3 x 3";
                end
             elseif(iscell(epsilon) && iscell(mu) && iscell(ksi) && iscell(zita))
                    res="";
             end
         catch
             try
                load(filename,'epsilon','mu','ksi','zita');
                medium=app.TModel.Domains(app.CurrentDomainIndex).Medium;
                if(ismatrix(epsilon) && ismatrix(mu) && ismatrix(ksi) && ismatrix(zita))
                    if(size(epsilon,1)==3 && size(epsilon,2)==3 && size(mu,1)==3 && size(mu,2)==3 && size(ksi,1)==3 && size(ksi,2)==3 && size(zita,1)==3 && size(zita,2)==3)
                        medium=medium.NoDispersive();medium=medium.Bianisotropic();
                        medium.Epsilon=epsilon;medium.Mu=mu;medium.Ksi=ksi;medium.Zita=zita;
                        app.TModel.Domains(app.CurrentDomainIndex).Medium=medium;
                        res = "Non Dispersive Bianisotropic Medium  for Domain "+ num2str(app.CurrentDomainIndex) +" Sucessfully loaded";
                    else
                        res=" Load Unsuccessful Tensor Matrices not 3 x 3";
                    end
                elseif(iscell(epsilon) && iscell(mu) && iscell(ksi) && iscell(zita))
                        res="";
                end
             catch
                 try
                     load(filename,'epsilon','mu','tag');
                    medium=app.TModel.Domains(app.CurrentDomainIndex).Medium;medium.Tag=tag;
                    if(isscalar(epsilon) && isscalar(mu))
                        medium=medium.NoDispersive();medium=medium.Isotropic();
                        medium.Epsilon=epsilon;medium.Mu=mu;
                        app.TModel.Domains(app.CurrentDomainIndex).Medium=medium;
                        res = "Non Dispersive Isotropic Medium  for Domain "+ num2str(app.CurrentDomainIndex) +" Sucessfully loaded ";
                    elseif(ismatrix(epsilon) && ismatrix(mu))
                        if(size(epsilon,1)==3 && size(epsilon,2)==3 && size(mu,1)==3 && size(mu,2)==3 )
                            medium=medium.NoDispersive();medium=medium.Anisotropic();
                            medium.Epsilon=epsilon;medium.Mu=mu;
                            app.TModel.Domains(app.CurrentDomainIndex).Medium=medium;
                            res = "Non Dispersive Anisotropic Medium  for Domain "+ num2str(app.CurrentDomainIndex) +" Sucessfully loaded ";
                        else
                            res=" Load Unsuccessful Tensor Matrices not 3 x 3";
                        end
                    elseif(iscell(epsilon) && iscell(mu))
                        res="";
                    else
                    end
                 catch
                     load(filename,'epsilon','mu');
                    medium=app.TModel.Domains(app.CurrentDomainIndex).Medium;
                    if(isscalar(epsilon) && isscalar(mu))
                        medium=medium.NoDispersive();medium=medium.Isotropic();
                        medium.Epsilon=epsilon;medium.Mu=mu;
                        app.TModel.Domains(app.CurrentDomainIndex).Medium=medium;
                        res = "Non Dispersive Isotropic Medium  for Domain "+ num2str(app.CurrentDomainIndex) +" Sucessfully loaded ";
                    elseif(ismatrix(epsilon) && ismatrix(mu))
                        if(size(epsilon,1)==3 && size(epsilon,2)==3 && size(mu,1)==3 && size(mu,2)==3 )
                            medium=medium.NoDispersive();medium=medium.Anisotropic();
                            medium.Epsilon=epsilon;medium.Mu=mu;
                            app.TModel.Domains(app.CurrentDomainIndex).Medium=medium;
                            res = "Non Dispersive Anisotropic Medium  for Domain "+ num2str(app.CurrentDomainIndex) +" Sucessfully loaded ";
                        else
                            res=" Load Unsuccessful Tensor Matrices not 3 x 3";
                        end
                    end
                 end

             end
         end

    end
end




