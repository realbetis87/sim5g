classdef Excitation
 %=========================================================================
 %{
    
    1. Type            : "POR","DIR"
    2. Edges           : Edge Structure updated for 2D Problem
    3. Facets          : Facet Structure updated for 2D Problem
    4. Vertices        : Vertex Structure updated for 2d Problem
    5. Freq            : frequency (Hz)
    6. MF              : MultiFrequency object
    7. Et              : Transverse Component Electric Field
                            vector (single frequency)
                            MF.NFx1 cell (MultiFrequency) 
    8. En              : Normal Component Electric Field
                            vector (single frequency)
                            MF.NFx1 cell (MultiFrequency) 
    9. Ht              : Transverse Component Magnetic Field 
                            vector (single frequency)
                            MF.NFx1 cell (MultiFrequency) 
    10.Hn              : Normal Component Magnetic Field
                            vector (single frequency)
                            MF.NFx1 cell (MultiFrequency) 
    11.Beta            : propagation constant (real part)
                            scalar (single frequency)
                            MFx1 vector (MultiFrequency)
    12.Alpha           : propagation constant (imaginary part)
                            scalar (single frequency)
                            MFx1 vector (MultiFrequency)
    13.LineBoundaries  : Indices to the Line Boundaries in the Excitation
                         2D Domain
    14.BoundaryIndices : Indices to the Surface Boundaries in the
                         Excitation 2D Domain
  
 %}
 %=========================================================================
    properties
        Type;
        Et;En;Ht;Hn;
        Alpha;Beta;Zita;
        MF=[];Freq=[];
        Assembly;Vector;KVector;
        Vertices=[];Edges=[];Facets=[];
        BoundaryIndices;LineBoundaries;
    end
    
    methods
        function obj = Excitation(varargin)
            if(nargin==0),return;
            elseif(nargin==1),obj.Type=varargin{1};
            elseif(nargin==2),obj.Type=varargin{1};
                if(isscalar(varargin{2})),obj.Freq=varargin{2};
                else,obj.MF=varargin{2};
                     obj.Et=cell(obj.MF.NF,1); obj.En=cell(obj.MF.NF,1);
                     obj.Ht=cell(obj.MF.NF,1); obj.Hn=cell(obj.MF.NF,1);
                     obj.Alpha=zeros(obj.MF.NF,1); obj.Beta=zeros(obj.MF.NF,1);
                end
            end
        end
        function obj = AddVertex(obj,VertexIndex),obj.Vertices(end+1)=VertexIndex;end
        function obj = AddEdge(obj,EdgeIndex),obj.Edges(end+1)=EdgeIndex;end
        function obj = AddFacet(obj,FacetIndex),obj.Facets(end+1)=FacetIndex;end
    end
end

