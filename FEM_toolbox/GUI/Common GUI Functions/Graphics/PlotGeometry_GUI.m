function [] = PlotGeometry_GUI(app,Mode)
        Geometry=app.TModel.model.Geometry;Vertices=Geometry.Vertices;xs=Vertices(:,1);ys=Vertices(:,2);zs=Vertices(:,3);
        xmin=min(xs);xmax=max(xs);ymin=min(ys);ymax=max(ys);zmin=min(zs);zmax=max(zs);xD=xmax-xmin;yD=ymax-ymin;zD=zmax-zmin;
        if(Mode==0),pdegplot(Geometry,"CellLabels","on","FaceLabels","on","FaceAlpha",0.5);title("")
        elseif(Mode==1),pdegplot(Geometry,"CellLabels","on","FaceAlpha",0.5);title("")
        elseif(Mode==2),pdegplot(Geometry,"FaceLabels","on","FaceAlpha",0.5);title("")
        elseif(Mode==3),pdegplot(Geometry,"FaceAlpha",0.5);title("")
        end
        hFig=gcf;hAx=gca;set(hAx,'YColor',[1.00,0.90,0.80]);set(hAx,'ZColor',[1.00,0.90,0.80]);set(hAx,'XColor',[1.00,0.90,0.80]);
        copyobj(allchild(hAx),app.UIAxes);view(app.UIAxes,15,25);   
        if(app.equalCheckBox.Value),axis(app.UIAxes,'equal');end
        if(app.tightCheckBox.Value),axis(app.UIAxes,'tight');end,close(hFig);
        app.UIAxes.XLim=[xmin-xD/10 xmax+xD/10];app.UIAxes.YLim=[ymin-yD/10 ymax+yD/10];app.UIAxes.ZLim=[zmin-zD/10 zmax+zD/10];
        set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
end

