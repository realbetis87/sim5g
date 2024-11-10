classdef Facet
%==========================================================================
%{
    Properties:
    
    1. Index         : Index of Facet
    2. Vertices            : Vector of Facet Vector Indices (3x1)
    3. OnBoundary          : Index of Boundary
    4. InElement           : Vector of Element Indices
    5. Id                  : 0 - Uknown
                             1 -  
    6. UknownIndex         : Index of DoF
    7. KnownIndex          : Known Numbering Index (Dirichlet Excitation)
    8. OpposingVertex      : (scalar or 2x1 vector) of opposing vertex in
                             tetrahedral element
    9. OnPort              : (index of Port)
    10.Surface             : Facet's Area
    11.SecondaryIndex      :  
    12.Barycenter          : [x,y,z] Barycentric of Facet
    13.Support             : Hierarchical Matrix Support
    14.Index2D             : 2D Problem Index
    15.Vertices2D          : (3x1) vector of Facet's 2D Vertex Indices
    16.Edges               : Indices of Facet Edges (used in 2D Problem)
    17.EdgeSigns           : Edge Signs (used in 2D Problem)
    18.Medium2D            : Medium (used in 2D Problem)
    19.PPair               : Index of corresponding facet in Periodic Pair  
    20.NormalVector        : Normal Vector to the facet's Surface.

%}
%==========================================================================
    properties
       Index;Vertices=zeros(3,1);OnBoundary=[];InElement=[];Id=0;KnownIndex=[];UknownIndex=[];OpposingVertex=[];OnPort=false;Surface;Barycenter=[];SecondaryIndex;Support;
       Index2D;Edges=zeros(3,1);EdgeSigns=zeros(3,1);Vertices2D;Medium2D;PPair;NormalVector;
    end
    
    methods
        function obj = Facet(varargin)
            if(nargin==0),return;
            elseif(nargin==6),Structure=varargin{1};obj.Index=varargin{2};obj.Vertices(1)=varargin{3};obj.Vertices(2)=varargin{4};obj.Vertices(3)=varargin{5};obj.InElement(1)=varargin{6};obj=obj.UpdateFacet(Structure);
            end
        end
        function obj = UpdateFacet(obj,Structure),vertex1=Structure.Vertices(obj.Vertices(1));vertex2=Structure.Vertices(obj.Vertices(2));vertex3=Structure.Vertices(obj.Vertices(3));
            obj.Surface=FacetSurface(vertex1,vertex2,vertex3);obj.Barycenter=FacetBarycenter(vertex1,vertex2,vertex3);
            obj.OpposingVertex(1)=setdiff([Structure.Elements(obj.InElement).Vertices],obj.Vertices);
        end
        function Structure = UpdateStructureFacet(obj,Structure),vertex1=Structure.Vertices(obj.Vertices(1));vertex2=Structure.Vertices(obj.Vertices(2));vertex3=Structure.Vertices(obj.Vertices(3));
            %-------------------------  Find Adjacent Element--------------
            CommonElements=intersect(vertex1.InElement,vertex2.InElement);CommonElements=intersect(CommonElements,vertex3.InElement);
            for ii=1:numel(CommonElements)
                if(CommonElements(ii)~=obj.InElement(1)),obj.InElement(2)=CommonElements(ii);
                    [sign,position] = FacetPositionSign(Structure.Elements(CommonElements(ii)),vertex1,vertex2,vertex3);
                    Structure.Elements(CommonElements(ii))=Structure.Elements(CommonElements(ii)).AddFacet(position,obj.Index,sign);
                    obj.OpposingVertex(2)=setdiff([Structure.Elements(CommonElements(ii)).Vertices],obj.Vertices);
                end
            end
            %------------------------------ Find Boundary -----------------

            Structure.Facets(obj.Index)=obj;
        end 
        function Structure = UpdateFacetBoundary(obj,Structure),edge1=Structure.Edges(obj.Edges(1));edge2=Structure.Edges(obj.Edges(2));edge3=Structure.Edges(obj.Edges(3));
                        B1=edge1.OnBoundary;B2=edge2.OnBoundary;B3=edge3.OnBoundary;BC1=intersect(B1,B2);BC=intersect(BC1,B3);
            if(numel(BC)>1)
                pause(0.1);
            end
            if(BC~=0),obj.OnBoundary=BC;Structure.Boundaries(BC)=Structure.Boundaries(BC).AddFacet(obj.Index);
            end
        end
    end
end
%----------------------- Local Functions ----------------------------------
function  A = FacetSurface(Vertex1,Vertex2,Vertex3),L1=EdgeLength(Vertex1,Vertex2);L2=EdgeLength(Vertex1,Vertex3);L3=EdgeLength(Vertex2,Vertex3);SemiPerimeter=(L1+L2+L3)/2;A=sqrt(SemiPerimeter*(SemiPerimeter-L1)*(SemiPerimeter-L2)*(SemiPerimeter-L3));end
function  L = EdgeLength(Vertex1,Vertex2),L=sqrt((Vertex1.X - Vertex2.X)^2 + (Vertex1.Y - Vertex2.Y)^2 + (Vertex1.Z - Vertex2.Z)^2);end
function  B = FacetBarycenter(Vertex1,Vertex2,Vertex3), B=[(Vertex1.X+Vertex2.X+Vertex3.X)/3;(Vertex1.Y+Vertex2.Y+Vertex3.Y)/3;(Vertex1.Z+Vertex2.Z+Vertex3.Z)/3;];end
function [Sign,Position] = FacetPositionSign(obj,Vertex1,Vertex2,Vertex3),Sign=0;Position=1;V=[Vertex1.Index Vertex2.Index Vertex3.Index];FacetPositions=[1 2 3 4;2 4 4 2;3 3 1 1;];ElementVertices=obj.Vertices;
    for ii=1:4,ElementV=ElementVertices(FacetPositions(:,ii));Check=double(eq(ElementV,V));if(nnz(Check)==3),Position=ii;Permutations=nnz(eye(3)-Check);if(nnz(Permutations)==0 || nnz(Permutations)==6),Sign=1;else,Sign=-1;end,end,end
end


%{
classdef Facet
%==========================================================================
%{
                            Properties
    1.Index
    2.Vertices
    3.Id
    4.UknownIndex
    5.KnownIndex
    6.OnBoundary
    7.InElement
    8.OpposingVertex
    9.Surface
    10.OnPort
    11.Barycenter
    12.SecondaryIndex

                              Constructor
    Facet(Index,Vertices)
%}
%==========================================================================
    properties
        Index;Vertices;Id=0;OnBoundary=0;KnownIndex=0;UknownIndex=0;InElement=[];OpposingVertex=[];OnPort=false;Surface;Barycenter=[];SecondaryIndex;
    end
    methods
        function obj = Facet(varargin)
            if(nargin==4),obj.Index=varargin{1};obj.Vertices(1)=varargin{2};obj.Vertices(2)=varargin{3};obj.Vertices(3)=varargin{4};
             obj=CheckBoundaryFacet(obj);obj=FacetNumbering(obj);obj=CheckElements(obj);
            end
        end
        function obj = CheckBoundaryFacet(obj),global BoundaryConditions;BoundaryFound=false;BoundaryCount=numel(BoundaryConditions);ii=0;
            while(BoundaryFound==false),ii=ii+1;boundary=BoundaryConditions(ii);if(OnGeometry(boundary.Geometry,obj)),obj.Id=boundary.Id;obj.OnBoundary=ii;BoundaryConditions(ii)=AddObjectOnBoundary(boundary,obj);BoundaryFound=true;end,if(ii+1>BoundaryCount),BoundaryFound=true;end,end
        end
        function obj = CheckElements(obj),global Elements;global Vertices;vertices=Vertices(obj.Vertices);CommonElements=intersect(vertices(1).InElement,vertices(2).InElement);CommonElements=intersect(CommonElements,vertices(3).InElement);
             for ii=1:numel(CommonElements),element=Elements(CommonElements(ii));[sign,pos]=FacetPositionSign(element,vertices(1),vertices(2),vertices(3));Elements(CommonElements(ii))=AddObjectOnElement(element,obj,sign,pos);obj=AddElementOnObject(obj,CommonElements(ii));end
             obj.Surface=FacetSurface(vertices(1),vertices(2),vertices(3));obj.Barycenter=[(vertices(1).X+vertices(2).X+vertices(3).X)/3;(vertices(1).Y+vertices(2).Y+vertices(3).Y)/3;(vertices(1).Z+vertices(2).Z+vertices(3).Z)/3;];
        end
        function obj = FacetNumbering(obj),global BoundaryConditions;global KnownCounter;global UknownCounter;
            if(obj.Id==0),UknownCounter=UknownCounter+1;obj.UknownIndex=UknownCounter;
            elseif(BoundaryConditions(obj.OnBoundary).Type=="Dirichlet"),KnownCounter=KnownCounter+1;obj.KnownIndex=KnownCounter;
           % elseif(obj.Id=='g' || obj.Id=='e' || obj.Id=='a'),UknownCounter=UknownCounter+1;obj.UknownIndex=UknownCounter;
            elseif(obj.Id=='g' || obj.Id=='a'),UknownCounter=UknownCounter+1;obj.UknownIndex=UknownCounter;
            elseif(obj.Id=='p'),UknownCounter=UknownCounter+1;obj.UknownIndex=UknownCounter;obj.OnPort=true;
            end
        end
    end
end
function  A = FacetSurface(Vertex1,Vertex2,Vertex3),L1=EdgeLength(Vertex1,Vertex2);L2=EdgeLength(Vertex1,Vertex3);L3=EdgeLength(Vertex2,Vertex3);SemiPerimeter=(L1+L2+L3)/2;A=sqrt(SemiPerimeter*(SemiPerimeter-L1)*(SemiPerimeter-L2)*(SemiPerimeter-L3));end
function  L = EdgeLength(Vertex1,Vertex2),L=sqrt((Vertex1.X - Vertex2.X)^2 + (Vertex1.Y - Vertex2.Y)^2 + (Vertex1.Z - Vertex2.Z)^2);end
function [Sign,Position] = FacetPositionSign(obj,Vertex1,Vertex2,Vertex3),Sign=0;Position=1;V=[Vertex1.Index Vertex2.Index Vertex3.Index];FacetPositions=[1 2 3 4;2 4 4 2;3 3 1 1;];ElementVertices=obj.Vertices;
    for ii=1:4,ElementV=ElementVertices(FacetPositions(:,ii));Check=double(eq(ElementV,V));if(nnz(Check)==3),Position=ii;Permutations=nnz(eye(3)-Check);if(nnz(Permutations)==0 || nnz(Permutations)==6),Sign=1;else,Sign=-1;end,end,end
end

%}