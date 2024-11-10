function [] = Plot2DDomain_GUI(app),boundary1=app.mainApp.TModel.Boundaries(app.BoundaryIndices(1));
    switch abs(boundary1.Axis)
        case 1
            for ii=1:numel(app.BoundaryIndices),boundary=app.mainApp.TModel.Boundaries(app.BoundaryIndices(ii));
                for jj=1:numel(boundary.Lines),lineBoundary=app.mainApp.TModel.LineBoundaries(boundary.Lines(jj));
                    for kk=1:numel(lineBoundary.Edges),edge=app.mainApp.TModel.Edges(lineBoundary.Edges(kk));
                        vertex1=app.mainApp.TModel.Vertices(edge.Vertices(1));vertex2=app.mainApp.TModel.Vertices(edge.Vertices(2));
                        if(app.LineBoundaryIndices(app.selectedLine)==lineBoundary.Index),plot([vertex1.Y vertex2.Y],[vertex1.Z vertex2.Z],'r','LineWidth',2);hold on;else,plot([vertex1.Y vertex2.Y],[vertex1.Z vertex2.Z],'k','LineWidth',2);hold on;end
                    end
                end
            end,xlabel("y axis");ylabel("z axis");
            hFig=gcf;hAx=gca;set(hAx,'YColor','r');set(hAx,'ZColor','r');set(app.UIAxes,'XColor','r');copyobj(allchild(hAx),app.UIAxes);close(hFig);
            set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
            app.UIAxes.XLabel.String="Y";app.UIAxes.YLabel.String="Z";
            axis(app.UIAxes,'equal');axis(app.UIAxes,'padded');
        case 2
            for ii=1:numel(app.BoundaryIndices),boundary=app.mainApp.TModel.Boundaries(app.BoundaryIndices(ii));
                for jj=1:numel(boundary.Lines),lineBoundary=app.mainApp.TModel.LineBoundaries(boundary.Lines(jj));
                    for kk=1:numel(lineBoundary.Edges),edge=app.mainApp.TModel.Edges(lineBoundary.Edges(kk));
                        vertex1=app.mainApp.TModel.Vertices(edge.Vertices(1));vertex2=app.mainApp.TModel.Vertices(edge.Vertices(2));
                        if(app.LineBoundaryIndices(app.selectedLine)==lineBoundary.Index),plot([vertex1.X vertex2.X],[vertex1.Z vertex2.Z],'r','LineWidth',2);hold on;else,plot([vertex1.X vertex2.X],[vertex1.Z vertex2.Z],'k','LineWidth',2);hold on;end
                    end
                end
            end,xlabel("x axis");ylabel("z axis");
            hFig=gcf;hAx=gca;set(hAx,'YColor','r');set(hAx,'ZColor','r');set(app.UIAxes,'XColor','r');copyobj(allchild(hAx),app.UIAxes);close(hFig);
            set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
            app.UIAxes.XLabel.String="X";app.UIAxes.YLabel.String="Z";
            axis(app.UIAxes,'equal');axis(app.UIAxes,'padded');
        case 3
            for ii=1:numel(app.BoundaryIndices),boundary=app.mainApp.TModel.Boundaries(app.BoundaryIndices(ii));
                for jj=1:numel(boundary.Lines),lineBoundary=app.mainApp.TModel.LineBoundaries(boundary.Lines(jj));
                    for kk=1:numel(lineBoundary.Edges),edge=app.mainApp.TModel.Edges(lineBoundary.Edges(kk));
                        vertex1=app.mainApp.TModel.Vertices(edge.Vertices(1));vertex2=app.mainApp.TModel.Vertices(edge.Vertices(2));
                        if(app.LineBoundaryIndices(app.selectedLine)==lineBoundary.Index),plot([vertex1.X vertex2.X],[vertex1.Y vertex2.Y],'r','LineWidth',2);hold on;else,plot([vertex1.X vertex2.X],[vertex1.Y vertex2.Y],'k','LineWidth',2);hold on;end
                    end
                end
            end,xlabel("x axis");ylabel("y axis");
            hFig=gcf;hAx=gca;set(hAx,'YColor','r');set(hAx,'ZColor','r');set(app.UIAxes,'XColor','r');copyobj(allchild(hAx),app.UIAxes);close(hFig);
            set(app.UIAxes,'YColor',[1.00,0.90,0.80]);set(app.UIAxes,'ZColor',[1.00,0.90,0.80]);set(app.UIAxes,'XColor',[1.00,0.90,0.80]);
            app.UIAxes.XLabel.String="X";app.UIAxes.YLabel.String="Y";
            axis(app.UIAxes,'equal');axis(app.UIAxes,'padded');
    end
end



%{

case 1
            %for ii=1:numel(app.BoundaryIndices),boundary=app.mainApp.TModel.Boundaries(app.BoundaryIndices(ii));
            %    for jj=1:numel(boundary.Edges),edge=app.mainApp.TModel.Edges(boundary.Edges(jj));
            %        vertex1=app.mainApp.TModel.Vertices(edge.Vertices(1));vertex2=app.mainApp.TModel.Vertices(edge.Vertices(2));
            %        plot([vertex1.Y vertex2.Y],[vertex1.Z vertex2.Z],'b');hold on;
            %    end
            %end
case 2
            %for ii=1:numel(app.BoundaryIndices),boundary=app.mainApp.TModel.Boundaries(app.BoundaryIndices(ii));
            %    for jj=1:numel(boundary.Edges),edge=app.mainApp.TModel.Edges(boundary.Edges(jj));
            %        vertex1=app.mainApp.TModel.Vertices(edge.Vertices(1));vertex2=app.mainApp.TModel.Vertices(edge.Vertices(2));
            %        plot([vertex1.X vertex2.X],[vertex1.Z vertex2.Z],'b');hold on;
            %    end
            %end
case 3
            %for ii=1:numel(app.BoundaryIndices),boundary=app.mainApp.TModel.Boundaries(app.BoundaryIndices(ii));
            %    for jj=1:numel(boundary.Edges),edge=app.mainApp.TModel.Edges(boundary.Edges(jj));
            %        vertex1=app.mainApp.TModel.Vertices(edge.Vertices(1));vertex2=app.mainApp.TModel.Vertices(edge.Vertices(2));
            %        plot([vertex1.X vertex2.X],[vertex1.Y vertex2.Y],'b');hold on;
            %    end
            %end
%}