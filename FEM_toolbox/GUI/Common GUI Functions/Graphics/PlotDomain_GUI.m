function [] = PlotDomain_GUI(app,domainIndex),cla(app.UIAxes,'reset');PlotGeometry_GUI(app,1);
%{
        cla(app.UIAxes,'reset');TModel=app.TModel;domain=TModel.Domains(domainIndex);
        Geometry=TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
        xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
        pdegplot(Geometry,"CellLabels","on","FaceAlpha",0.5);hold on;
        for ii=1:numel(domain.Elements),element=TModel.Elements(domain.Elements(ii));
            v1=TModel.Vertices(element.Vertices(1));v2=TModel.Vertices(element.Vertices(2));v3=TModel.Vertices(element.Vertices(3));v4=TModel.Vertices(element.Vertices(4));
            P1=[v1.X v2.X v3.X v4.X];P2=[v1.Y v2.Y v3.Y v4.Y];P3=[v1.Z v2.Z v3.Z v4.Z];fill3(P1,P2,P3,[0.4660 0.6740 0.1880],'FaceAlpha',0.6,'LineStyle','none');hold on;
            P4=[v4.X v3.X v2.X v1.X];P5=[v4.Y v3.Y v2.Y v1.Y];P6=[v4.Z v3.Z v2.Z v1.Z];fill3(P4,P5,P6,[0.4660 0.6740 0.1880],'FaceAlpha',0.6,'LineStyle','none');hold on;
        end
        hFig=gcf;hAx=gca;set(hAx,'YColor',[1.00,0.90,0.80]);set(hAx,'ZColor',[1.00,0.90,0.80]);set(hAx,'XColor',[1.00,0.90,0.80]);
        copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);
        if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
        if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
        app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];
        set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
%}
end

