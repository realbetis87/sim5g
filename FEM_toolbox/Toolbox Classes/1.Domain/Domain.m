classdef Domain
 %=========================================================================
 %{
    
    Properties:

        1. Index                : Domain Index    
        2. Dim                  : 2 (2D Domain), 3 (3D Domain) 
        3. Medium               : Medium Object
        5. Elements             : Vector of Element Indices in Domain
        6. Vertices             : Vector of Vertex Indices in Domain
        7. ExternalBoundaries   : Vector of Boundary Indices on the exterior of the Domain
        8. InternalBoundaries   : Vector of Boundary Indices on the interior of the Domain
    
    Methods:
           
        1. obj = Domain()                         : Empty Object Constructor
        2. obj = Domain(Domain Index)     
        3. obj = Domain(Domain Index,ExternalBoundaries,InternalBoundaries)

        
 %}
 %======================================================================
    properties
        Index;Dim=3;Medium=Medium(0);Elements=[];Vertices=[];ExternalBoundaries;InternalBoundaries;IsEmpty=false;
    end
    
    methods
        function obj = Domain(varargin)
            if(nargin==0),return;
            elseif(nargin==1),obj.Index=varargin{1};
            elseif(nargin==3),obj.Index=varargin{1};obj.ExternalBoundaries=varargin{2};obj.InternalBoundaries=varargin{3};
            end
        end
        function obj = AddElement(obj,ElementIndex),obj.Elements(end+1)=ElementIndex;end
        function obj = AddVertex(obj,VertexIndex),obj.Vertices(end+1)=VertexIndex;end
        function obj = DefineMedium(obj,Material),obj.Medium=Material;end
      end
end

