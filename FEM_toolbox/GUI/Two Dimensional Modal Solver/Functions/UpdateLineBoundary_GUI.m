function [] = UpdateLineBoundary_GUI(app),Plot2DDomain_GUI(app);
    LineBoundary=app.TModel.LineBoundaries(app.LineBoundaryIndices(app.selectedLine));
    switch LineBoundary.Type
        case "PEC",app.BoundaryTypeSelection.ValueIndex=2;app.BoundaryTypeButton.Enable=false;app.BoundaryTypeButton.Visible=false;
        case "CON",app.BoundaryTypeSelection.ValueIndex=1;app.BoundaryTypeButton.Enable=false;app.BoundaryTypeButton.Visible=false;
        case "PMC",app.BoundaryTypeSelection.ValueIndex=3;app.BoundaryTypeButton.Enable=false;app.BoundaryTypeButton.Visible=false;
        case "GRA",app.BoundaryTypeSelection.ValueIndex=5;app.BoundaryTypeButton.Enable=true;app.BoundaryTypeButton.Visible=true;app.BoundaryTypeButton.Text="σ";
        case "ABC",app.BoundaryTypeSelection.ValueIndex=4;app.BoundaryTypeButton.Enable=false;app.BoundaryTypeButton.Visible=false;
    end
end

