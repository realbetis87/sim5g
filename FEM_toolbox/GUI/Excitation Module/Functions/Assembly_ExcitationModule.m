function [] = Assembly_ExcitationModule(app)
    %--------------- Assembly Initialization ------------------------------
    if(app.TModel.Frequency.NF==1),Assembly=FEMAssembly("Excitation","3D");Assembly.Bian=BianCheck(app.TModel);Assembly.Dispersive=false;
        Assembly.NE=ReturnNE(app.TModel);Assembly.N=ReturnNB(app.TModel);Assembly.NB=Assembly.N-Assembly.NE;
        Assembly.IsDir=DirCheck(app.TModel);
        if(Assembly.IsDir), Assembly.KNE=ReturnKNE(app.TModel);Assembly.KN=ReturnKNB(app.TModel);Assembly.KNB=Assembly.KN-Assembly.KNE;end
        Assembly.IsPort=PortCheck(app.TModel);
        app.TModel.Assembled=Assembly;
    else,Assembly=FEMAssembly("Excitation","3D",app.TModel.Frequency.NF);Assembly.Dispersive=true;
        if(BianCheck(app.TModel)),Assembly.Bian=true;Assembly=Assembly.Init3DExcitationBian_MF();else,Assembly=Assembly.Init3DExcitationMode_MF();end
        app.TModel.Assembled=Assembly;
    end

    %-------------- Scaling Selection -------------------------------------
    switch(app.ScalingSelection.ValueIndex)
        case 1
        case 2
        case 3
        case 4
        case 5
        case 6
    end
    app.TModel=Excitation_Assembly(app.TModel);
end

function [BianFlag] = BianCheck(TModel),BianFlag=false;for ii=1:numel(TModel.Domains),dom=TModel.Domains(ii);med=dom.Medium;if(med.Type=="Bian"),BianFlag=true;end,end,end
function [DirCheck] = DirCheck(TModel),DirCheck=false;for ii=1:numel(TModel.Boundaries),bou=TModel.Boundaries(ii);if(bou.Type=="DIR"),DirCheck=true;end,end,end
function [PCheck] = PortCheck(TModel),PCheck=false;for ii=1:numel(TModel.Boundaries),bou=TModel.Boundaries(ii);if(bou.Type=="POR"),PCheck=true;end,end,end
function [NE] = ReturnNE(TModel),edges=[TModel.Edges];NE=max([edges.UknownIndex]);end
function [NB] = ReturnNB(TModel),facets=[TModel.Facets];NB=max([facets.UknownIndex]);end
function [KNE] = ReturnKNE(TModel),edges=[TModel.Edges];KNE=max([edges.KnownIndex]);end
function [KNB] = ReturnKNB(TModel),facets=[TModel.Facets];KNB=max([facets.KnownIndex]);end