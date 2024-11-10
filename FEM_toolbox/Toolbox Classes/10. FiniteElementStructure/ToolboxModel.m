classdef ToolboxModel
%==========================================================================
%{
        Container Class for all Finite Element Structures 
    1.  model                :   Matlab PDEmodel object (contains Geometry & Mesh)
    2.  Domains              :   Vector of Domain objects
    3.  Boundaries           :   Vector of Boundary (Surface) objects
    4.  LineBoundaries       :   Vector of Boundary (Line) objects
    5.  Vertices             :   Vector of Vertex Objects
    6.  Elements             :   Vector of Element Objects
    7.  Edges                :   Vector of Edge Objects
    8.  Facets               :   Vector of Facet Objects
    9.  Media                :   Vector of Medium Objects
    10. Assembled            :   FEMAssembly Object
    11. Assembled_2D         :   FEMAssembly Object
    11. Solution             :   FEMSolution Object
    12. Cond                 :   0 - Nothing Done,1 - Geometry Done,2 - FEM Initializations Done,3 - Frequency Done
                                 4 - Domains Done,5 - Boundaries Done,6 - Assembly Done,7 - Solution Done
    13.Boundary_Excitations  :  Vector of Excitation Objects
    14.UknownExcVector       
    15.KnownExcVector
%}
%==========================================================================
    properties
        model=[];Domains=[];Boundaries=[];Elements=[];Vertices=[];Edges=[];Facets=[];LineBoundaries=[];Media=[];Cond=0;
        Frequency;Domain2D;
        NumberOfDomains=0;NumberOfBoundaries=0;NumberOfLineBoundaries=0;NumberOfElements=0;NumberOfVertices=0;NumberOfEdges=0;NumberOfFacets=0;NumberOfBoundaryVertices=0;
        Assembled=[];Assembled_2D=[];Solution=[];Boundary_Excitations=[];
        UknownExcVector=[];KnownExcVector=[];
        
    end
    
    methods
        function obj = ToolboxModel(varargin)
            if(nargin==0),return;
            elseif(nargin==4),obj.model=varargin{1};obj.Boundaries=varargin{2};obj.LineBoundaries=varargin{3};obj.Domains=varargin{4};
                              obj.NumberOfBoundaries=size(obj.Boundaries,2);obj.NumberOfDomains=size(obj.Domains,2);obj.NumberOfLineBoundaries=size(obj.LineBoundaries,2);obj.Media=Medium.empty(0,obj.NumberOfDomains);
            elseif(nargin==5),obj.Domains=varargin{1};obj.LineBoundaries=varargin{2};obj.Vertices=varargin{3};obj.Edges=varargin{4};obj.Facets=varargin{5};obj.Media=Medium.empty(0,obj.NumberOfDomains);
            elseif(nargin==8)
                obj.model=varargin{1};obj.Domains=varargin{2};obj.Boundaries=varargin{3};
                obj.Elements=varargin{4};obj.Vertices=varargin{5};obj.Edges=varargin{6};obj.Facets=varargin{7};
                obj.NumberOfBoundaryVertices=varargin{8};
                obj.NumberOfDomains=numel(obj.Domains);obj.NumberOfBoundaries=numel(obj.Boundaries);
                obj.NumberOfElements=numel(obj.Elements);obj.NumberOfVertices=numel(obj.Vertices);
                obj.NumberOfEdges=numel(obj.Edges);obj.NumberOfFacets=numel(obj.Facets);obj.Media=Medium.empty(0,obj.NumberOfDomains);
            end
        end
        function obj = AddFrequency(obj,Frequency),obj.Frequency=Frequency;end      
        function [obj,Index] = AddExcitation(obj,ExcitationObject)
            if(isempty(obj.Boundary_Excitations)),Index=1;obj.Boundary_Excitations=ExcitationObject;
            else,temp=obj.Boundary_Excitations;Ne=numel(temp);obj.Boundary_Excitations=Excitation.empty(Ne+1,0);
                for ii=1:Ne,obj.Boundary_Excitations(ii)=temp(ii);end,obj.Boundary_Excitations(Ne+1)=ExcitationObject;Index=Ne+1;
            end
        end
        function obj = RemoveExcitation(obj,ExcitationIndex)
            Ne=numel(obj.Boundary_Excitations);temp=obj.Boundary_Excitations;obj.Boundary_Excitations=Excitation.empty(Ne-1,0);
            for ii=1:ExcitationIndex-1,obj.Boundary_Excitations(ii)=temp(ii);end
            for jj=ExcitationIndex+1:Ne,obj.Boundary_Excitations(jj-1)=temp(jj);end
        end
    end
end

