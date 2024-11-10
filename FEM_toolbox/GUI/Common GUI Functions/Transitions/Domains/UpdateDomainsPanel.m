function [] = UpdateDomainsPanel(app),PlotGeometry_GUI(app,1),app.DomainSelection.Items={};
    for ii=1:app.TModel.NumberOfDomains,str="Domain "+int2str(ii); app.DomainSelection.Items{end+1}=char(str);end,app.DomainFlags=zeros(app.TModel.NumberOfDomains,1);
    pause(0.1)
end

