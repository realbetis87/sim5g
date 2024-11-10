function [] = TensorBoundary_GUI(app)
    if(app.TensorCheckBox.Value),app.TModel.Boundaries(app.CurrentBoundaryIndex).Tensorial();
    else,app.TModel.Boundaries(app.CurrentBoundaryIndex).Scalar();
    end
end

