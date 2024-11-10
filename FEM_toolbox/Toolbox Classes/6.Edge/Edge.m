classdef Edge
%==========================================================================
%{
    Properties

        1. Index               : Index of Edge
        2. UknownIndex               : Index of DoF 
        3. KnownIndex                : Known Numbering Index (Dirichlet Excitation)
        4. SecondaryIndex            : (Redundant)
        5. Vertices                  : Vector of Edge Vector Indices (2x1)
        6. InElement                 : Vector Of Edge's Element Indices 
        7. OnBoundary                : Index of Boundary
        8. OnLine                    : Index of Line Boundary
        9. Id                        : if on Boundary ~ inherits Boundary Id
                                       else 0 
        10.PPair                     : Index of Periodic Pair
        11.Length                    : Edge Length
        12.Support                   : Hierarchical Matrices Support
        13.MidPoint                  : Hierarchical Matrices Midpoint
        14.Index2D                   : Index of Edge in 2D Problem
        15.IndexE                    : Index of Electric Field DoF in 2D Problem
        16.IndexH                    : Index of Magnetic Field DoF in 2D Problem
        17.InFacet                   : Vector of Edge's Facet 2D Indices
        18.Vertices2D                : (2 x 1) vector of Edge's Vertices in
                                        2D

    Constructor
%}
%==========================================================================
    
    properties
        Vertices;Index;UknownIndex;KnownIndex;SecondaryIndex;InElement=[];OnBoundary=[];OnLine=[];Id;Length;Support;MidPoint;
        Index2D;IndexE;IndexH;InFacet;Vertices2D;PPair;
    end
    
    methods
        function obj = Edge(varargin)
            if(nargin==0),return;
            elseif(nargin==1),obj.Index=varargin{1};
            elseif(nargin==5),obj.Index=varargin{2};obj.Vertices(1)=varargin{3};obj.Vertices(2)=varargin{4};obj.InElement(1)=varargin{5};
                Structure=varargin{1};obj=obj.LengthAndBoundary(Structure);
            end
        end
        function obj = AddInElement(obj,inElement),obj.InElement(end+1)=inElement;end
        function obj = LengthAndBoundary(obj,Structure),Vertex1=Structure.Vertices(obj.Vertices(1));Vertex2=Structure.Vertices(obj.Vertices(2));
                 obj.Length=sqrt((Vertex1.X - Vertex2.X)^2 + (Vertex1.Y - Vertex2.Y)^2 + (Vertex1.Z-Vertex2.Z)^2);obj.MidPoint=[(Vertex1.X+Vertex2.X)/2;(Vertex1.Y+Vertex2.Y)/2;(Vertex1.Z+Vertex2.Z)/2;];
                 LV1=Vertex1.OnLine;LV2=Vertex2.OnLine;BV1=Vertex1.OnBoundary;BV2=Vertex2.OnBoundary;LV=intersect(LV1,LV2);BV=intersect(BV1,BV2);
                 obj.OnLine=LV;
                 if(isscalar(BV)),obj.OnBoundary=BV;
                 elseif(~isempty(BV))
                     if(~isempty(obj.OnLine)),obj.OnBoundary=BV;
                     else,NF=nearestFace(Structure.model.Geometry,obj.MidPoint');if(isscalar(NF)),obj.OnBoundary=NF;else,BV=intersect(BV,NF);obj.OnBoundary=BV;end
                     end
                 end
        end
    end
end

