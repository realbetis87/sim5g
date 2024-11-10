function [] = InitializationsComplete_EigenModeModule(app)
    %---------------- Left Panel ------------------------------------------
    PopulateBoundariesTree(app);PopulateDomainsTree(app);app.LoadGeometryButton.Text="Continue";
    app.LoadGrid.RowHeight={'0x','1x','0x'};app.ImportGeometryGrid.RowHeight={'1x','3x','3x'};
    %---------------- Central Panel ---------------------------------------
    app.CentralGrid.RowHeight={'1x','0x','0x','10x','1x'};
    app.BoundariesTreePanel.Visible=true;app.BoundariesTreePanel.Enable=true;
    app.DomainsTreePanel.Visible=true;app.DomainsTreePanel.Enable=true;
    UpdateDomainsPanel(app);app.GeometryButton.FontColor=[0.47,0.67,0.19];app.GeometryButton.Enable=true;
    %---------------- Plot Menu -------------------------------------------
    app.PlotButton.Enable=true;app.PlotSelectionArea.Enable=true;app.PlotSelection.Enable=true;   
    app.PlotSelectionArea.Items{1}='Plot All Domains';app.PlotSelectionArea.Items{2}='Plot All Boundaries';
    app.PlotSelection.Enable=true;app.PlotSelection.Visible=true;
    for ii=1:app.TModel.model.Geometry.NumCells,str="Plot Domain : " + int2str(ii);app.PlotSelectionArea.Items{end+1}=char(str);end
    for ii=1:app.TModel.model.Geometry.NumFaces,str="Plot Boundary : "+int2str(ii);app.PlotSelectionArea.Items{end+1}=char(str);end
    app.SaveButton.Enable=true;app.SaveButton.Visible=true;
end

function [] = PopulateBoundariesTree(app)
    p=uitreenode(app.DomainsTree,"Text","Domains");
    for ii=1:app.TModel.model.Geometry.NumCells,domain=app.TModel.Domains(ii);
        pp=uitreenode(p,"Text","SubDomain "+int2str(ii));
        pp1=uitreenode(pp,"Text","Number Of Elements "+int2str(numel(domain.Elements)));
        pp2=uitreenode(pp,"Text","Number Of Vertices "+int2str(numel(domain.Vertices)));
    end
end
function [] = PopulateDomainsTree(app)
    q=uitreenode(app.BoundariesTree,"Text","Boundaries");
    for jj=1:app.TModel.model.Geometry.NumFaces,boundary=app.TModel.Boundaries(jj);
       qq=uitreenode(q,"Text","Boundary "+int2str(jj));
       qq2=uitreenode(qq,"Text","Number Of Vertices "+int2str(numel(boundary.Vertices)));
       qq3=uitreenode(qq,"Text","Number Of Edges "+int2str(numel(boundary.Edges)));
       qq4=uitreenode(qq,"Text","Number Of Facets "+int2str(numel(boundary.Facets)));
    end
end