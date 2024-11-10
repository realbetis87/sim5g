function [] = TwoDModalSolver_StartUp(app)
    [app.TModel,app.BoundaryIndices,app.Excitation_Res,app.LineBoundaryIndices]=Prepare2DBoundaries(app.TModel,app.mainApp.CurrentBoundaryIndex);app.Equation=app.Excitation_Res.Assembly.Equation;
    %----------------------------------------------------------------------
    Plot2DDomain_GUI(app);
    %----------------------------------------------------------------------
    for ii=1:numel(app.LineBoundaryIndices),app.BoundarySelection.Items{end+1}=char(int2str(ii));end,app.selectedLine=1;
    app.ResultsPanel.Enable=false;app.ResultsPanel.Visible=false;
    %----------------------------------------------------------------------
    if(app.mainApp.TModel.Boundaries(app.mainApp.CurrentBoundaryIndex).Type=="POR"),app.Type="POR";else,app.Type="DIR";end,app.Excitation_Res.Type=app.Type;
    NewMessage(app,"2D Modal Solver initilized for " + num2str(numel(app.BoundaryIndices)) + "Boundaries and " + num2str(numel(app.LineBoundaryIndices)) + " line boundaries");
    NewMessage(app,"2D Mesh consists of " + num2str(numel(app.Excitation_Res.Vertices)) + " Vertices, " +  num2str(numel(app.Excitation_Res.Edges)) + " Edges and " +  num2str(numel(app.Excitation_Res.Facets)) +" Facets.");  
    switch app.Equation
        case 1,NewMessage(app,"Equations Selected : E-H (Maxwell)Tangential Normal - Tangential Normal");app.EquationLabel_1.Text="E-H Formulation";app.EquationLabel_2.Text="2D 1/2";
        case 2,NewMessage(app,"Equations Selected : E-H (Maxwell) Tangential - Tangential");app.EquationLabel_1.Text="E-H Formulation";app.EquationLabel_2.Text="2D";
        case 3,NewMessage(app,"Equation Selected :  E - (Vector Wave Equation) Tangential Normal ");app.EquationLabel_1.Text="E Formulation";app.EquationLabel_2.Text="2D 1/2";
    end
    if(app.TModel.Frequency.NF~=1),app.eigenvalues_selected=zeros(app.TModel.Frequency.NF,1);end
    UpdateLineBoundary_GUI(app);PrepareTwoDAssembly(app);
end


function [] = PrepareTwoDAssembly(app)
    if(app.TModel.Frequency.NF>1)
        app.Excitation_Res.Assembly.Matrix_A=cell(app.TModel.Frequency.NF,1);
        app.Excitation_Res.Assembly.Matrix_B=cell(app.TModel.Frequency.NF,1);
    end
end