function [toolboxModel] = DetermineBoundaryGeometry(toolboxModel),error=1000*eps;
    for ib=1:numel(toolboxModel.Boundaries);boundary=toolboxModel.Boundaries(ib);
        boundaryVertices=[toolboxModel.Vertices(boundary.Vertices)];Xs=[boundaryVertices.X];Ys=[boundaryVertices.Y];Zs=[boundaryVertices.Z];
        if(all(abs(Xs-Xs(1))<=error)),boundary.Position=Xs(1);bBox=toolboxModel.model.Geometry.boundingBox();xmin=bBox(1,1);xmax=bBox(1,2);
            if(abs(Xs(1)-xmin)<=error),boundary.Axis=1;boundary=boundary.PEC();
            elseif(abs(Xs(1)-xmax)<=error),boundary.Axis=-1;boundary=boundary.PEC();
            else,boundary.Axis=1.5;
            end
        elseif(all(abs(Ys-Ys(1))<=error)),boundary.Position=Ys(1);bBox=toolboxModel.model.Geometry.boundingBox();ymin=bBox(2,1);ymax=bBox(2,2);
            if(abs(Ys(1)-ymin)<=error),boundary.Axis=2;boundary=boundary.PEC();
            elseif(abs(Ys(1)-ymax)<=error),boundary.Axis=-2;boundary=boundary.PEC();
            else,boundary.Axis=2.5;
            end
        elseif(all(abs(Zs-Zs(1))<=error)),boundary.Position=Zs(1);bBox=toolboxModel.model.Geometry.boundingBox();zmin=bBox(3,1);zmax=bBox(3,2);
            if(abs(Zs(1)-zmin)<=error),boundary.Axis=3;boundary=boundary.PEC();
            elseif(abs(Zs(1)-zmax)<=error),boundary.Axis=-3;boundary=boundary.PEC();
            else,boundary.Axis=3.5;
            end
        else,boundary.Axis=4;
        end,toolboxModel.Boundaries(ib)=boundary;
    end
end

