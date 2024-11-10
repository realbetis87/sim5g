function [] = PlotFEMGeometryMesh(app)
       if(app.PlotSelection.Value=="Vertices")
            if(app.PlotSelectionArea.ValueIndex==1),PlotAllVertices(app);
            elseif(app.PlotSelectionArea.ValueIndex==2),PlotAllBoundaryVertices(app);
            elseif(app.PlotSelectionArea.ValueIndex>2 && app.PlotSelectionArea.ValueIndex<=2+app.TModel.NumberOfDomains),domain=app.PlotSelectionArea.ValueIndex-2;PlotVerticesOnDomain(app,domain);
            elseif(app.PlotSelectionArea.ValueIndex>2+app.TModel.NumberOfDomains),boundary=app.PlotSelectionArea.ValueIndex-2-app.TModel.NumberOfDomains;PlotVerticesOnBoundary(app,boundary);
            end
       elseif(app.PlotSelection.Value=="Edges")
            if(app.PlotSelectionArea.ValueIndex==1),PlotAllEdges(app);
            elseif(app.PlotSelectionArea.ValueIndex==2),PlotAllBoundaryEdges(app);   
            elseif(app.PlotSelectionArea.ValueIndex>2 && app.PlotSelectionArea.ValueIndex<=2+app.TModel.NumberOfDomains),domain=app.PlotSelectionArea.ValueIndex-2;PlotEdgesOnDomain(app,domain);
            elseif(app.PlotSelectionArea.ValueIndex>2+app.TModel.NumberOfDomains),boundary=app.PlotSelectionArea.ValueIndex-app.TModel.NumberOfDomains-2;
                PlotEdgesOnBoundary(app,boundary);
            end
        elseif(app.PlotSelection.Value=="Facets")
            if(app.PlotSelectionArea.ValueIndex==1),PlotAllFacets(app);
            elseif(app.PlotSelectionArea.ValueIndex==2)   
            elseif(app.PlotSelectionArea.ValueIndex>2 && app.PlotSelectionArea.ValueIndex<=2+app.TModel.NumberOfDomains),domain=app.PlotSelectionArea.ValueIndex-2;PlotFacetsOnDomain(app,domain);
            elseif(app.PlotSelectionArea.ValueIndex>2+app.TModel.NumberOfDomains),boundary=app.PlotSelectionArea.ValueIndex-app.TModel.NumberOfDomains-2;
                PlotFacetsOnBoundary(app,boundary);
            end
        end
end
%------------------------------ Plot Vertices -----------------------------
function [] = PlotVerticesOnBoundary(app,boundary),cla(app.UIAxes,'reset');
    pdegplot(app.TModel.model.Geometry,"CellLabels","on","FaceLabels","on","FaceAlpha",0.5);hold on;
    for ii=1:numel(app.TModel.Vertices),if(app.TModel.Vertices(ii).OnBoundary==boundary),plot3(app.TModel.Vertices(ii).X,app.TModel.Vertices(ii).Y,app.TModel.Vertices(ii).Z,'o','Color','k');hold on;end,end
    hFig=gcf;hAx=gca;copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);
   if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
   if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
   set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
   Geometry=app.TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
   xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
   app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];
end
function [] = PlotAllBoundaryVertices(app)
   pdegplot(app.TModel.model.Geometry,"CellLabels","on","FaceLabels","on","FaceAlpha",0.5);hold on;
   for ii=1:numel(app.TModel.Vertices),if(~isempty(app.TModel.Vertices(ii).OnBoundary)),plot3(app.TModel.Vertices(ii).X,app.TModel.Vertices(ii).Y,app.TModel.Vertices(ii).Z,'o','Color','k');hold on;end,end
   hFig=gcf;hAx=gca;copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);
   if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
   if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
   set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
   Geometry=app.TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
   xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
   app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];

end
function [] = PlotVerticesOnDomain(app,domain),cla(app.UIAxes,'reset');TModel=app.TModel;dom=TModel.Domains(domain);
    pdegplot(app.TModel.model.Geometry,"CellLabels","on","FaceLabels","on","FaceAlpha",0.5);hold on;
    for ii=1:numel(dom.Vertices),vert=TModel.Vertices(dom.Vertices(ii));
        plot3(vert.X,vert.Y,vert.Z,'o','Color','k','MarkerFaceColor','k');hold on;
    end
   hFig=gcf;hAx=gca;copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);
   if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
   if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
   set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
   Geometry=app.TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
   xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
   app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];
end
function [] = PlotAllVertices(app),cla(app.UIAxes,'reset');
    pdegplot(app.TModel.model.Geometry,"CellLabels","on","FaceLabels","on","FaceAlpha",0.5);hold on;
    for ii=1:numel(app.TModel.Vertices),plot3(app.TModel.Vertices(ii).X,app.TModel.Vertices(ii).Y,app.TModel.Vertices(ii).Z,'o','Color','k');hold on;end
    hFig=gcf;hAx=gca;copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);
   if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
   if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
   set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
   Geometry=app.TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
   xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
   app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];
end
%----------------------------- Plot Edges ---------------------------------
function [] = PlotAllEdges(app),cla(app.UIAxes,'reset');
    pdegplot(app.TModel.model.Geometry,"CellLabels","on","FaceLabels","on","FaceAlpha",0.5);hold on;
    for ii=1:numel(app.TModel.Edges),Vertex1=app.TModel.Vertices(app.TModel.Edges(ii).Vertices(1));Vertex2=app.Vertices(app.TModel.Edges(ii).Vertices(2));plot3([Vertex1.X Vertex2.X],[Vertex1.Y Vertex2.Y],[Vertex1.Z Vertex2.Z],'Color','k');hold on;end
    hFig=gcf;hAx=gca;copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);
   if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
   if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
   set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
   Geometry=app.TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
   xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
   app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];
end
function [] = PlotEdgesOnDomain(app,domain),cla(app.UIAxes,'reset');TModel=app.TModel;dom=TModel.Domains(domain);
    pdegplot(app.TModel.model.Geometry,"CellLabels","on","FaceLabels","on","FaceAlpha",0.5);hold on;
    for ii=1:numel(dom.Elements),el=TModel.Elements(dom.Elements(ii));
        for jj=1:6,edge=TModel.Edges(el.Edges(jj));
            Vertex1=TModel.Vertices(edge.Vertices(1));Vertex2=TModel.Vertices(edge.Vertices(2));
            plot3([Vertex1.X Vertex2.X],[Vertex1.Y Vertex2.Y],[Vertex1.Z Vertex2.Z],'Color','k');hold on;
        end
    end
   hFig=gcf;hAx=gca;copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);   
   if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
   if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
   set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
   Geometry=app.TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
   xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
   app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];
end
function [] = PlotEdgesOnBoundary(app,boundary),cla(app.UIAxes,'reset');
    Edges=app.TModel.Edges;pdegplot(app.TModel.model.Geometry,"CellLabels","on","FaceLabels","on","FaceAlpha",0.5);hold on;
    for ii=1:numel(Edges)
        if(Edges(ii).OnBoundary==boundary),Vertex1=app.TModel.Vertices(Edges(ii).Vertices(1));Vertex2=app.TModel.Vertices(Edges(ii).Vertices(2));
            plot3([Vertex1.X Vertex2.X],[Vertex1.Y Vertex2.Y],[Vertex1.Z Vertex2.Z],'Color','k');hold on;
        end
    end
    hFig=gcf;hAx=gca;copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);    
   if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
   if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
   set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
   Geometry=app.TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
   xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
   app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];
    
end
%----------------------------- Plot Facets --------------------------------
function [] = PlotAllFacets(app),cla(app.UIAxes,'reset');
    pdegplot(app.TModel.model.Geometry,"CellLabels","on","FaceLabels","on","FaceAlpha",0.5);hold on;
    for ii=1:numel(app.TModel.Facets),Vertex1=app.TModel.Vertices(app.Facets(ii).Vertices(1));Vertex2=app.TModel.Vertices(app.Facets(ii).Vertices(2));Vertex3=app.TModel.Vertices(app.Facets(ii).Vertices(3));
        P1=[Vertex1.X,Vertex2.X,Vertex3.X];P2=[Vertex1.Y,Vertex2.Y,Vertex3.Y];P3=[Vertex1.Z,Vertex2.Z,Vertex3.Z];fill3(P1,P2,P3,[0.8500 0.3250 0.0980],'FaceAlpha',0.4);hold on;
    end
   hFig=gcf;hAx=gca;copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);
   if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
   if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
   set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
   Geometry=app.TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
   xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
   app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];
end
function [] = PlotFacetsOnBoundary(app,boundary),cla(app.UIAxes,'reset');
    pdegplot(app.TModel.model.Geometry,"CellLabels","on","FaceLabels","on","FaceAlpha",0.5);hold on;
    for ii=1:numel(app.TModel.Facets)
        if(app.TModel.Facets(ii).OnBoundary==boundary)
            Vertex1=app.TModel.Vertices(app.TModel.Facets(ii).Vertices(1));Vertex2=app.TModel.Vertices(app.TModel.Facets(ii).Vertices(2));Vertex3=app.TModel.Vertices(app.TModel.Facets(ii).Vertices(3));
            P1=[Vertex1.X,Vertex2.X,Vertex3.X];P2=[Vertex1.Y,Vertex2.Y,Vertex3.Y];P3=[Vertex1.Z,Vertex2.Z,Vertex3.Z];fill3(P1,P2,P3,[0.8500 0.3250 0.0980],'FaceAlpha',0.4);hold on;
        end
    end
   hFig=gcf;hAx=gca;copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);
   if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
   if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
   set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
   Geometry=app.TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
   xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
   app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];
end
function [] = PlotFacetsOnDomain(app,domain),cla(app.UIAxes,'reset');TModel=app.TModel;dom=TModel.Domains(domain);
    pdegplot(app.TModel.model.Geometry,"CellLabels","on","FaceLabels","on","FaceAlpha",0.5);hold on;
    for ii=1:numel(dom.Elements),el=TModel.Elements(dom.Elements(ii));
        for jj=1:4,facet=TModel.Facets(el.Facets(jj));
            Vertex1=TModel.Vertices(facet.Vertices(1));Vertex2=TModel.Vertices(facet.Vertices(2));Vertex3=TModel.Vertices(facet.Vertices(3));
            P1=[Vertex1.X,Vertex2.X,Vertex3.X];P2=[Vertex1.Y,Vertex2.Y,Vertex3.Y];P3=[Vertex1.Z,Vertex2.Z,Vertex3.Z];fill3(P1,P2,P3,[0.8500 0.3250 0.0980],'FaceAlpha',0.4);hold on;
        end
    end
   hFig=gcf;hAx=gca;copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);   
   if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
   if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
   set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
   Geometry=app.TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
   xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
   app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];
end