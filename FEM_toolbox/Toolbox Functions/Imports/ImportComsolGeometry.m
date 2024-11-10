%--------------------------------------------------------------------------
%{
    Toolbox function importing Comsol Mesh
    Comsol Mesh Data : p        : matrix (Number of Mesh Vertices x 3)
                       t        : matrix (Number of Mesh Tetrahedral Elements x 4)
                       meshdata : comsol structure 
    1. Load Mesh
    2. Initialize Matlab PDE 3D Model from Mesh
    3. Find External Surface and Line Boundaries
    4. Initialize LineBoundary, Boundary and Domain Structures
%}
%--------------------------------------------------------------------------
function [TModel,res] = ImportComsolGeometry(filename)
    try 
        %---------------------- Load Meshdata and Mesh --------------------
        load(filename,'p','t','meshdata');t=double(t);elem=meshdata.elementity{2};model=createpde(3);
        try
             Geometry = geometryFromMesh(model,p,t,elem);res=0;
             %----------- Computational Domain Exterior Boundaries--------
             [ExteriorBoundaries,ExteriorEdges] = ReturnExternalBoundaries(filename);
             %------------------- Initialize Toolbox Model-----------------
             Domains=Domain.empty(0,Geometry.NumCells);Boundaries=Boundary.empty(0,Geometry.NumFaces);
             LineBoundaries=Boundary.empty(0,Geometry.NumEdges);
             TModel=ToolboxModel(model,Boundaries,LineBoundaries,Domains);   
             %------------------------ Line Boundaries --------------------
             for ii=1:Geometry.NumEdges,attachedFaces=facesAttachedToEdges(Geometry,ii);
                if(nnz(ismember(ExteriorEdges,ii))==1),TModel.LineBoundaries(ii)=Boundary(ii,1,true,attachedFaces);
                else,TModel.LineBoundaries(ii)=Boundary(ii,1,false,attachedFaces);
                end
            end
            %-------------------------- Boundaries ------------------------
            for ii=1:Geometry.NumFaces,FaceEdges=faceEdges(Geometry,ii);
                if(nnz(ismember(ExteriorBoundaries,ii))==1),TModel.Boundaries(ii)=Boundary(ii,2,true,FaceEdges);
                else,TModel.Boundaries(ii)=Boundary(ii,2,false,FaceEdges);
                end
            end
            %---------------------- Domains -------------------------------------------
            for ii=1:Geometry.NumCells,ExternalBoundaries=cellFaces(Geometry,ii,"external");InternalBoundaries=cellFaces(Geometry,ii,"internal");
                TModel.Domains(ii)=Domain(ii,ExternalBoundaries,InternalBoundaries);TModel.Domains(ii).Medium=Medium("Medium " + num2str(ii));
            end
        catch
            TModel=[];res=1;
        end
        
    catch
        TModel=[];res=2;
    end
end

function [ExteriorBoundaries,ExteriorEdges] = ReturnExternalBoundaries(filename)
    model=createpde(3);load(filename,'p','t','meshdata');t=double(t);elem=meshdata.elementity{2};
     MergedGeometry = geometryFromMesh(model,p,t,elem);
    while MergedGeometry.NumCells>1
        MergedGeometry=mergeCells(MergedGeometry,[1 2]);
    end
    ExteriorBoundaries=cellFaces(MergedGeometry,1,"external");
    ExteriorEdges=cellEdges(MergedGeometry,1,"external");
end