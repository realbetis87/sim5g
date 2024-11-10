function [] = Assembly_EigenModeModule(app)
    %--------------- Assembly Initialization ------------------------------
    if(app.TModel.Frequency.NF==1),Assembly=FEMAssembly("EigenMode","3D");Assembly.Bian=BianCheck(app.TModel);app.TModel.Assembled=Assembly;
    else,Assembly=FEMAssembly("EigenMode","3D",app.TModel.Frequency.NF);Assembly.Dispersive=true;Assembly.E_Scaling=1;Assembly.B_Scaling=1;
        if(BianCheck(app.TModel)),Assembly.Bian=true;Assembly=Assembly.Init3DEigenModeBian_MF();else,Assembly=Assembly.Init3DEigenMode_MF();end
         app.TModel.Assembled=Assembly;
    end
    if(app.NCheck.Value),app.TModel.Assembled.EigenValue="n";
    else,app.TModel.Assembled.EigenValue="k";
    end
    %----------------------------------------------------------------------
    edges=[app.TModel.Edges];facets=[app.TModel.Facets];
    NE=max([edges.UknownIndex]);N=max([facets.UknownIndex]);NB=N-NE;
    app.TModel.Assembled.NE=NE;app.TModel.Assembled.NB=NB;app.TModel.Assembled.N=N;
    app.NE_Info.Text=num2str(NE);app.NB_Info.Text=num2str(NB);app.N_Info.Text=num2str(N);
    %-------------- Scaling Selection -------------------------------------
    switch(app.ScalingSelection.ValueIndex)
        case 1,app.TModel.Assembled.E_Scaling=1;app.TModel.Assembled.B_Scaling=1;
        case 2,app.TModel.Assembled.E_Scaling=1;app.TModel.Assembled.B_Scaling=0;
        case 3,app.TModel.Assembled.E_Scaling=0;app.TModel.Assembled.B_Scaling=1;
        case 4,app.TModel.Assembled.E_Scaling=0;app.TModel.Assembled.B_Scaling=0;
    end
    %------------------- Propagation Axis ---------------------------------
    switch app.PropagationAxisSelection.ValueIndex
        case 1,Assembly=app.TModel.Assembled;Assembly.PropagationAxis="x";app.TModel.Assembled=Assembly;
        case 2,Assembly=app.TModel.Assembled;Assembly.PropagationAxis="y";app.TModel.Assembled=Assembly;
        case 3,Assembly=app.TModel.Assembled;Assembly.PropagationAxis="z";app.TModel.Assembled=Assembly;
    end
    %-------------------- Call Assembly -----------------------------------
    pause(0.1);app.TModel=EigenMode_Assembly(app.TModel);
end

function [BianFlag] = BianCheck(TModel),BianFlag=false;for ii=1:numel(TModel.Domains),dom=TModel.Domains(ii);med=dom.Medium;if(med.Type=="Bian"),BianFlag=true;end,end,end
