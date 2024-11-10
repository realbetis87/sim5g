function [] = UpdateBoundariesPanel(app)
    app.ColumnMainGrid.ColumnWidth={'1x','4x','1x'};BoundariesInfo=cell(app.TModel.NumberOfBoundaries,3);
    for ii=1:app.TModel.NumberOfBoundaries,boundary=app.TModel.Boundaries(ii);
        BoundariesInfo{ii,1}=char(num2str(ii));BoundariesInfo{ii,3}=char(boundary.Type);
        switch boundary.Axis
            case 1  ,BoundariesInfo{ii,2}='+x';
            case -1 ,BoundariesInfo{ii,2}='-x';
            case 1.5,BoundariesInfo{ii,2}='x';
            case 2  ,BoundariesInfo{ii,2}='+y';
            case -2 ,BoundariesInfo{ii,2}='-y';
            case 2.5,BoundariesInfo{ii,2}='+y';
            case 3,BoundariesInfo{ii,2}='+z';
            case -3,BoundariesInfo{ii,2}='-z';
            case 3.5,BoundariesInfo{ii,2}='z';
            case 4,BoundariesInfo{ii,2}='n';
        end
        str=int2str(ii);app.CurrentBoundarySelection.Items{end+1}=char(str);
    end,app.BoundaryFlags=zeros(app.TModel.NumberOfBoundaries,1);app.BoundariesInfo=BoundariesInfo;app.BoundariesTable.Data=BoundariesInfo;
    UpdateCurrentBoundary(app);
end