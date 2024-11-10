function [] = PDBoundarySelection_StartUp(app)
    boundary=app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex);ax=boundary.Axis;pos=boundary.Position;
    switch boundary.Axis
        case 1,if(sign(boundary.Axis)>0),app.OrLabel.Text="+x";else,app.OrLabel.Text="-x";end
        case 2,if(sign(boundary.Axis)>0),app.OrLabel.Text="+y";else,app.OrLabel.Text="-y";end
        case 3,if(sign(boundary.Axis)>0),app.OrLabel.Text="+z";else,app.OrLabel.Text="-z";end
        case 4,app.OrLabel.Text="undetermined";
    end
    if(~iemspty(boundary.Pos)),app.PosLabel.Text=num2str(boundary.Pos);else,app.PosLabel.Text="undetermined"
    switch (boundary.Type)
        case "POR"

        case "DIR"
    end
end

