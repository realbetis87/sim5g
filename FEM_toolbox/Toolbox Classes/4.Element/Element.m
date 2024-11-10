classdef Element
%==========================================================================
%{
        
        Properties:

        1. Index   :  Index of element.
        2. SubDomain     :  Index of SubDomain containing the Element
        3. Type          :  "Tri" - 2D Triangular Element
                            "Tet" - 3D Tetrahedral Element
        4. Vertices      :  Vector of Global Indices of Element Vertices
                                (3x1) - Triangular Element
                                (4x1) - Tetrahedral Element
        5. Edges         :  Vector of Global Indices of Element Edges
                                (3x1) - Triangular Element
                                (6x1) - Tetrahedral Element
        6. EdgeSigns     :  Signs of Element.Edges 
        7. Facets        :   Vector of Global Index of Element Facets
                                null - Triangular Element
                                (4x1)- Tetrahedral Element
        8.FacetSigns     : Signs of Element.Facets
        9. Barycenter    : Barycentric Point
        10.Volume        : Element Volume
        11.As            : Simplex Coordinate Coefficients
        12.Bs            : Simplex Coordinate Coefficients
        13.Cs            : Simplex Coordinate Coefficients
        14.Ds            : Simplex Coordinate Coefficients
        Functions
        
        i.   Element()
        ii.  Element(GlobalIndex,Vertices)
%}
%==========================================================================
    properties
        Index;SubDomain;Type;Vertices=zeros(4,1);Edges=zeros(6,1);Facets=zeros(4,1);EdgeSigns=zeros(6,1);FacetSigns=zeros(4,1);Volume=0;Barycenter=zeros(3,1);As=zeros(3,1);Bs=zeros(3,1);Cs=zeros(3,1);Ds=zeros(3,1);
    end
    
    methods
        function obj = Element(varargin)
            if(nargin==0),return;
            elseif(nargin==1),obj.Index=varargin{1};
            elseif(nargin==2),obj.Index=varargin{1};obj.Vertices=varargin{2};
            elseif(nargin==3),obj.Index=varargin{1};obj.Vertices=varargin{2};obj.SubDomain=varargin{3};
            elseif(nargin==5)
            end
        end
        function obj = UpdateElement(obj,Vertices),vertices=Vertices(obj.Vertices);
            x=[vertices.X];y=[vertices.Y];z=[vertices.Z];
            De=det([1 x(1) y(1) z(1);1 x(2) y(2) z(2);1 x(3) y(3) z(3);1 x(4) y(4) z(4);]);obj.Volume=(1/6)*abs(De);
            a(1)=det ([1 x(1) y(1) z(1); 0 x(2) y(2) z(2); 0 x(3) y(3) z(3); 0 x(4) y(4) z(4)]) / De;b(1)=det ([1 1 y(1) z(1); 1 0 y(2) z(2); 1 0 y(3) z(3); 1 0 y(4) z(4)]) / De;c(1)=det ([1 x(1) 1 z(1); 1 x(2) 0 z(2); 1 x(3) 0 z(3); 1 x(4) 0 z(4)]) / De;d(1)=det ([1 x(1) y(1) 1; 1 x(2) y(2) 0; 1 x(3) y(3) 0; 1 x(4) y(4) 0]) / De; 
            a(2)=det ([0 x(1) y(1) z(1); 1 x(2) y(2) z(2); 0 x(3) y(3) z(3); 0 x(4) y(4) z(4)]) / De;b(2)=det ([1 0 y(1) z(1); 1 1 y(2) z(2); 1 0 y(3) z(3); 1 0 y(4) z(4)]) / De;c(2)=det ([1 x(1) 0 z(1); 1 x(2) 1 z(2); 1 x(3) 0 z(3); 1 x(4) 0 z(4)]) / De;d(2)=det ([1 x(1) y(1) 0; 1 x(2) y(2) 1; 1 x(3) y(3) 0; 1 x(4) y(4) 0]) / De; 
            a(3)=det ([0 x(1) y(1) z(1); 0 x(2) y(2) z(2); 1 x(3) y(3) z(3); 0 x(4) y(4) z(4)]) / De;b(3)=det ([1 0 y(1) z(1); 1 0 y(2) z(2); 1 1 y(3) z(3); 1 0 y(4) z(4)]) / De;c(3)=det ([1 x(1) 0 z(1); 1 x(2) 0 z(2); 1 x(3) 1 z(3); 1 x(4) 0 z(4)]) / De;d(3)=det ([1 x(1) y(1) 0; 1 x(2) y(2) 0; 1 x(3) y(3) 1; 1 x(4) y(4) 0]) / De; 
            a(4)=det ([0 x(1) y(1) z(1); 0 x(2) y(2) z(2); 0 x(3) y(3) z(3); 1 x(4) y(4) z(4)]) / De;b(4)=det ([1 0 y(1) z(1); 1 0 y(2) z(2); 1 0 y(3) z(3); 1 1 y(4) z(4)]) / De;c(4)=det ([1 x(1) 0 z(1); 1 x(2) 0 z(2); 1 x(3) 0 z(3); 1 x(4) 1 z(4)]) / De;d(4)=det ([1 x(1) y(1) 0; 1 x(2) y(2) 0; 1 x(3) y(3) 0; 1 x(4) y(4) 1]) / De; 
            obj.As=a;obj.Bs=b;obj.Cs=c;obj.Ds=d;obj.Barycenter=[sum(x)/4;sum(y/4);sum(z/4);];
        end
        function obj = AddEdge(obj,index,Edge,EdgeSign),obj.Edges(index)=Edge;obj.EdgeSigns(index)=EdgeSign;end
        function obj = AddFacet(obj,index,Facet,FacetSign),obj.Facets(index)=Facet;obj.FacetSigns(index)=FacetSign;end
      end
end
