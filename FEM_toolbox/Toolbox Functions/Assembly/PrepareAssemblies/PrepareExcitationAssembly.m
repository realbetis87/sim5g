function [TModel] = PrepareExcitationAssembly(TModel)
    Assembly=FEMAssembly("Excitation","3D");Assembly.E_Scaling=0;Assembly.B_Scaling=0;
    bianFlag=false;
    for ii=1:numel(TModel.Domains),domain=TModel.Domains(ii);medium=domain.Medium;
        if(medium.Type=="Bian"),bianFlag=true;end
    end,Assembly.Bian=bianFlag;dirFlag=false;
    for ii=1:numel(TModel.Boundary_Excitations),excitation=TModel.Boundary_Excitations(ii);
        if(excitation.Type=="DIR"),dirFlag=true;end
    end,Assembly.IsDir=dirFlag;
    TModel.Assembled=Assembly;
end

