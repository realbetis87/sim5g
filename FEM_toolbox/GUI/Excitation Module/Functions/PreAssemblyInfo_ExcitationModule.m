function [] = PreAssemblyInfo_ExcitationModule(app),TModel=app.TModel;
    NE= max([TModel.Edges.UknownIndex]);N=max([TModel.Facets.UknownIndex]);NB=N-NE;
    KNE=max([TModel.Edges.KnownIndex]);KN=max([TModel.Facets.KnownIndex]);KNB=KN-KNE;
    NP=0;for ii=1:numel(TModel.Boundary_Excitations),if(TModel.Boundary_Excitations(ii).Type=="POR"),NP=NP+numel(TModel.Boundary_Excitations.Edges);end,end
    app.TModel.Assembled.NE=NE;app.TModel.Assembled.NB=NB;app.TModel.Assembled.N=N;
    app.TModel.Assembled.KNE=KNE;app.TModel.Assembled.KNB=KNB;app.TModel.Assembled.KN=KN;
    app.NE_Info.Text=num2str(NE);app.NB_Info.Text=num2str(NB);app.N_Info.Text=num2str(N);
    app.KNE_Info.Text=num2str(KNE);app.KNB_Info.Text=num2str(KNB);app.KN_Info.Text=num2str(KN);
    app.PN_Info.Text=num2str(NP);
end