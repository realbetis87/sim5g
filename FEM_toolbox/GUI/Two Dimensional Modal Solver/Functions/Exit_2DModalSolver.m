function [] = Exit_2DModalSolver(app)
        excitation=app.Excitation_Res;
            if(isempty(app.TModel.Boundary_Excitations))
                for ii=1:numel(app.BoundaryIndices),app.TModel.Boundaries(app.BoundaryIndices(ii)).ExcitationIndex=1;end
                app.TModel.Boundary_Excitations=excitation;
                if(app.TModel.Frequency.NF~=1)
                    for ii=1:numel(app.BoundaryIndices)
                        app.TModel.Boundaries(app.BoundaryIndices(ii)).Dispersive=true;
                        app.TModel.Boundaries(app.BoundaryIndices(ii)).Param=zeros(app.TModel.Frequency.NF,1);
                        for ff=1:app.TModel.Frequency.NF
                            app.TModel.Boundaries(app.BoundaryIndices(ii)).Param(ff)=excitation.Beta(ff);
                        end
                     end
                else
                    for ii=1:numel(app.BoundaryIndices)
                        app.TModel.Boundaries(app.BoundaryIndices(ii)).Param=excitation.Beta;
                    end
                end
            else,n=numel(app.TModel.Boundary_Excitations);
                for ii=1:numel(app.BoundaryIndices),app.TModel.Boundaries(app.BoundaryIndices(ii)).ExcitationIndex=n+1;end
                temp=app.TModel.Boundary_Excitations;app.TModel.Boundary_Excitations=Excitation.empty(n+1,0);
                for jj=1:n,app.TModel.Boundary_Excitations(jj)=temp(jj);end,app.TModel.Boundary_Excitations(n+1)=excitation;
            end

            app.mainApp.TModel=app.TModel;UpdateCurrentBoundary(app.mainApp);
end

