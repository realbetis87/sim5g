classdef Vertex 
%==========================================================================
%{
    Properties:

    1. Index             : Index of Vertex
    2. InDomain          : Index of Vertex Domain
    3. X,Y,Z             : Vertex Coordinates
    4. InElement         : Vector of Vertex Element Indices
    5. OnBoundary        : Boundary Index - if Vertex on Bounary
                                0        - Vertex not on Boundary
    6. OnLine            : Line Boundary Index - if Vertex on Line Boundary
                                0              - Vertex not on Line Boundary
    7. OnExterior        : true/false
    8. Index2D           : Index of Vertex in 2D Problem
    9. InFacet           : Array of Vertex Facet Indices (2D Problem)
    10.IndexE            : Index of Electric Field DoF in 2D Problem
    11.IndexH            : Index of Magnetic Field DoF in 2D Problem
%}
%==========================================================================
    properties
        Index;InDomain;X;Y;Z;InElement=[];
        OnBoundary=0;OnLine=0;OnExterior=false;
        Index2D;InFacet;IndexE;IndexH;
    end
    methods
        function obj = Vertex(varargin)
            if(nargin==0),return;
            elseif(nargin==5),obj.Index=varargin{1};obj.X=varargin{2};obj.Y=varargin{3};obj.Z=varargin{4};obj.InDomain=varargin{5};
            elseif(nargin==7),structure=varargin{1};obj.Index=varargin{2};obj.X=varargin{3};obj.Y=varargin{4};obj.Z=varargin{5};obj.InDomain=varargin{6};obj.OnBoundary=varargin{7};obj.OnExterior=obj.IsExternal(structure);
            elseif(nargin==8),structure=varargin{1};obj.Index=varargin{2};obj.X=varargin{3};obj.Y=varargin{4};obj.Z=varargin{5};obj.InDomain=varargin{6};obj.OnBoundary=varargin{7};obj.OnLine=varargin{8};obj.OnExterior=obj.IsExternal(structure);
            end
        end
        function res = IsExternal(obj,structure),res=false;
            if(~isempty(obj.OnLine) & obj.OnLine~=0),res=structure.LineBoundaries(obj.OnLine(1)).OnExterior;
            elseif(~isempty(obj.OnBoundary) & obj.OnBoundary~=0),for ii=1:numel(obj.OnBoundary),if(structure.Boundaries(obj.OnBoundary(ii)).OnExterior),res=true;end,end
            end
        end
    end
end
