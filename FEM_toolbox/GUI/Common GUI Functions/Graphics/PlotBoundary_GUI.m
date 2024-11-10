function [] = PlotBoundary_GUI(app,boundaryIndex),cla(app.UIAxes,'reset');
    if(isvector(boundaryIndex) && numel(boundaryIndex)>1),TModel=app.TModel;
        Geometry=TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
        xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
        pdegplot(Geometry,"FaceLabels","on","FaceAlpha",0.5);hold on;
        for kk=1:numel(boundaryIndex),boundary=app.TModel.Boundaries(boundaryIndex(kk));
            for ii=1:numel(boundary.Facets),facet=app.TModel.Facets(boundary.Facets(ii));
                Vertex1=app.TModel.Vertices(facet.Vertices(1));Vertex2=app.TModel.Vertices(facet.Vertices(2));Vertex3=app.TModel.Vertices(facet.Vertices(3));
                P1=[Vertex1.X,Vertex2.X,Vertex3.X];P2=[Vertex1.Y,Vertex2.Y,Vertex3.Y];P3=[Vertex1.Z,Vertex2.Z,Vertex3.Z];fill3(P1,P2,P3,[0.8500 0.3250 0.0980],'FaceAlpha',0.6,'LineStyle','none');hold on;
            end
        end
        hFig=gcf;hAx=gca;set(hAx,'YColor',[1.00,0.90,0.80]);set(hAx,'ZColor',[1.00,0.90,0.80]);set(hAx,'XColor',[1.00,0.90,0.80]);
        copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);
        if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
        if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
        app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];
        set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
    else,boundary=app.TModel.Boundaries(boundaryIndex);TModel=app.TModel;
        Geometry=TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
        xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
        pdegplot(Geometry,"FaceLabels","on","FaceAlpha",0.5);hold on;
        if(boundary.Type=="PBC")
            for ii=1:numel(boundary.Facets),facet=app.TModel.Facets(boundary.Facets(ii));
                Vertex1=app.TModel.Vertices(facet.Vertices(1));Vertex2=app.TModel.Vertices(facet.Vertices(2));Vertex3=app.TModel.Vertices(facet.Vertices(3));
                P1=[Vertex1.X,Vertex2.X,Vertex3.X];P2=[Vertex1.Y,Vertex2.Y,Vertex3.Y];P3=[Vertex1.Z,Vertex2.Z,Vertex3.Z];fill3(P1,P2,P3,[0.8500 0.3250 0.0980],'FaceAlpha',0.6,'LineStyle','none');hold on;
            end
            if(~isempty(boundary.Param)),slave=app.TModel.Boundaries(boundary.Param);
                for ii=1:numel(slave.Facets),facet=app.TModel.Facets(slave.Facets(ii));
                    Vertex1=app.TModel.Vertices(facet.Vertices(1));Vertex2=app.TModel.Vertices(facet.Vertices(2));Vertex3=app.TModel.Vertices(facet.Vertices(3));
                    P1=[Vertex1.X,Vertex2.X,Vertex3.X];P2=[Vertex1.Y,Vertex2.Y,Vertex3.Y];P3=[Vertex1.Z,Vertex2.Z,Vertex3.Z];fill3(P1,P2,P3,[0.4660 0.6740 0.1880],'FaceAlpha',0.6,'LineStyle','none');hold on;
                end
            end
        else
            for ii=1:numel(boundary.Facets),facet=app.TModel.Facets(boundary.Facets(ii));
                Vertex1=app.TModel.Vertices(facet.Vertices(1));Vertex2=app.TModel.Vertices(facet.Vertices(2));Vertex3=app.TModel.Vertices(facet.Vertices(3));
                P1=[Vertex1.X,Vertex2.X,Vertex3.X];P2=[Vertex1.Y,Vertex2.Y,Vertex3.Y];P3=[Vertex1.Z,Vertex2.Z,Vertex3.Z];fill3(P1,P2,P3,[0.8500 0.3250 0.0980],'FaceAlpha',0.6,'LineStyle','none');hold on;
            end
        end
        hFig=gcf;hAx=gca;set(hAx,'YColor',[1.00,0.90,0.80]);set(hAx,'ZColor',[1.00,0.90,0.80]);set(hAx,'XColor',[1.00,0.90,0.80]);
        copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);
        if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
        if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
        app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];
        set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
    end
end

    
