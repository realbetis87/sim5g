function [] = PrepareExcitationAssembly_ExcitationModule(app)
    if(isempty(app.TModel.Assembled)),app.TModel.Assembled=FEMAssembly("Excitation","3D");end
    if(BianCheck(app.TModel)),app.TModel.Assembled.Bian=true;else,app.TModel.Assembled.Bian=false;end
    if(DirCheck(app.TModel)),app.TModel.Assembled.IsDir=true;else,app.TModel.Assembled.IsDir=false;end
    if(app.TModel.Frequency.NF>1),app.TModel.Assembled.Dispersive=true;app.TModel.Assembled.NF=app.TModel.Frequency.NF;app.TModel.Assembled=app.TModel.Assembled.Init3DExcitation_MF();
    else,app.TModel.Assembled.Dispersive=false;
    end
    if(PortCheck(app.TModel)),app.TModel.Assembled.IsPort=true;else,app.TModel.Assembled.IsPort=false;end
    switch app.ScalingSelection.ValueIndex
        case 1,app.TModel.Assembled.E_Scaling=1;app.TModel.Assembled.B_Scaling=1;
        case 2,app.TModel.Assembled.E_Scaling=1;app.TModel.Assembled.B_Scaling=0;
        case 3,app.TModel.Assembled.E_Scaling=0;app.TModel.Assembled.B_Scaling=1;
        case 4,app.TModel.Assembled.E_Scaling=0;app.TModel.Assembled.B_Scaling=0;
    end
    
end

function flag = BianCheck(TModel),flag=false;
    for ii=1:numel(TModel.Domains),medium=TModel.Domains(ii).Medium;
        if(medium.Type=="Bian"),flag=true;end
    end
end
function flag = PortCheck(TModel),flag=false;
       for ii=1:numel(TModel.Boundary_Excitations),excitation=TModel.Boundary_Excitations(ii);
            if(excitation.Type=="POR"),flag=true;end
       end
end
function flag = DirCheck(TModel),flag=false;
       for ii=1:numel(TModel.Boundary_Excitations),excitation=TModel.Boundary_Excitations(ii);
            if(excitation.Type=="DIR"),flag=true;end
       end
end
