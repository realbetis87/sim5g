classdef Boundary
%==========================================================================
%{
        Properties:
    
            A.Geometry Properties
    
        1.  Index                        : Boundary Index
        2.  Dim                          : 1 - Line Boundary
                                           2 - Surface Boundary
        3.  Axis                         : +\-1 - x Axis normal to the Boundary
        Boundary                            1.5 - x Axis normal to the Boundary, Boundary Internal
                                           +\-2 - y Axis normal to the Boundary 
                                            2.5 - y Axis normal to the Boundary, Boundary Internal
                                           +\-3 - z Axis normal to the Boundary
                                            3.5 - z Axis normal to the Boundary, Boundary Internal
                                            4   - not a plane boundary   
        4.  Position                     : x Coordinate if Boundary.Axis == 1
                                           y Coordinate if Boundary.Axis == 2    
                                           z Coordinate if Boundary.Axis == 3
        5.  Lines                        : Discrete Geometry Edges (indices) attached to Boundary
                                           (applicable for Surface Boundaries)
        6.  Faces                        : Discrete Geometry Faces (indices) attached to Boundary
                                           (applicable for Line Boundaries)
        7.  Exterior                     : true\false - Boundary is on the exterior of the computational domain.
  

            B. Mesh Properties
        
        7.  Vertices                     : Vector of Vertex Indices on Boundary
        8.  Edges                        : Vector of Edge Indices on Boundary 
        9.  Facets                       : Vector of Facet Indices on Boundary

            C. Boundary Condition Properites

        10. Type                         : "PEC"   - Perfect Electric Conductor   
                                           "PMC"   - Perfect Magnetic Conductor
                                           "DIR"   - Dirichlet Boundary Condition (Excitation)
                                           "POR"   - Port Boundary Condition (Excitation)
                                           "GRA"   - Graphene
                                           "IBC"   - Impedance Boundary Condition
                                           "PBC"   - Periodic Boundary Condition
                                           "ABB"   - Absorbing Boundary Condition with propagation constant (β) as parameter (Port without Excitation)
                                           "ABZ"   - Absorbing Boundary Condition with wave impedance (Z) as parameter (Port without Excitation)
                                           "CON"   - Continuity
                                           "ABC"   - Absorbing Boundary Condition
        11. Param                         :"ABB"   - propagation constant b
                                           "ABZ"   - wave impedance Z
                                           "GRA"   - graphene conductivity  
                                           "IBC"   - conductivity
                                           "POR"   - wave impedance
                                           "PBC"   - Boundary Pair
        12. Dispersive                    : true\false (dispersive Boundary Parameter)
        13. Master                        : true\false (PBC Boundary)
        14. Tensor                        : true\false (Applicable for ABZ boundary - wave impedance is tensorial)
        15. Id                            : Vertex,Edge,Facet Id
                                            0\1 (Uknown DoF, Known DoF)
                                              2 (PEC Id for EH 2D Formulations)
                                              3 (PMC Id for EH 2D Formulations)
        16.ExcitationIndex                : If Boundary Port\Dirichlet
                                            index to the
                                            TModel.Excitations()
        17.PortParamType                  : 0 Param == beta
                                            1 Param == Wave Impedance scalar
                                            2 Param == Wave Impedance Tensor
        18. PortType                      : 0  2D Modal Excitation
                                            1  Plane Wave Excitation
        19. PlaneWave                     : [Ex, Ey, Ez];
 %}
%==========================================================================
    properties
        Index;Dim=2;Axis=[];Position=[];Lines=[];Faces=[];OnExterior;
        Vertices=[];Edges=[];Facets=[];ExcitationIndex=[];
        Type="CON";Id=0;Param;Dispersive=false;Tensor=false;Master;Vec;
        PlaneWave=[0;0;0;];PortParamType=0;PortType=0;
    end
    methods
        function obj = Boundary(varargin)
            if(nargin==0)
            elseif(nargin==4),obj.Index=varargin{1};obj.Dim=varargin{2};obj.OnExterior=varargin{3};
                if(obj.Dim==1),obj.Faces=varargin{4};
                elseif(obj.Dim==2),obj.Lines=varargin{4};
                end
            end
        end
        function obj = AddEdge(obj,edge),obj.Edges(end+1)=edge;end
        function obj = AddFacet(obj,facet),obj.Facets(end+1)=facet;end
        function obj = PEC(obj),obj.Dispersive=false;obj.Tensor=false;obj.Param=[];obj.Type="PEC";obj.Id=1;end
        function obj = CON(obj),obj.Dispersive=false;obj.Tensor=false;obj.Param=[];obj.Type="CON";obj.Id=0;end
        function obj = PMC(obj),obj.Dispersive=false;obj.Tensor=false;obj.Param=[];obj.Type="PMC";obj.Id=0;end
        function obj = Graphene(varargin),obj=varargin{1};obj=obj.PEC();obj.Id=0;
            if(nargin==1),obj.Type="GRA";obj.Param=0;
            elseif(nargin==2),obj.Type="GRA";param=varargin{2};
                if(isa(param,"MultiFrequency")),obj.Param=zeros(param.NF,1);obj.Dispersive=true;
                elseif(isvector(param)),obj.Param=param;obj.Dispersive=true;
                else,obj.Param=param;
                end
            end
        end
        function obj = ABC_Z(varargin),obj=varargin{1};obj=obj.PEC();obj.Type="ABZ";obj.Id=0;
            if(nargin==1),obj.Param=0;
            elseif(nargin==2),param=varargin{2};obj.Type="ABZ";obj.Id=0;
                if(isa(param,"MultiFrequency")),obj.Param=zeros(param.NF,1);obj.Dispersive=true;
                elseif(isvector(param)),obj.Param=param;obj.Dispersive=true;
                else,obj.Param=param;
                end
            end
        end
        function obj = ABC_B(varargin),obj=varargin{1};obj=obj.PEC();obj.Type="ABB";obj.Id=0;
                if(nargin==1),obj.Param=0;
                elseif(nargin==2),param=varargin{2};obj.Type="ABB";obj.Id=0;
                    if(isa(param,"MultiFrequency")),obj.Param=zeros(param.NF,1);obj.Dispersive=true;
                    elseif(isvector(param)),obj.Param=param;obj.Dispersive=true;
                    else,obj.Param=param;
                    end
                end
        end
        function obj = IBC(varargin),obj=varargin{1};obj=obj.PEC();obj.Type="IBC";obj.Id=0;
            if(nargin==1),obj.Param=0;
            elseif(nargin==2),param=varargin{2};obj.Type="IBC";obj.Id=0;
                if(isa(param,"MultiFrequency")),obj.Param=zeros(param.NF,1);obj.Dispersive=true;
                elseif(isvector(param)),obj.Param=param;obj.Dispersive=true;
                else,obj.Param=param;
                end
            end
        end
        function obj = ABC_ZT(varargin),obj=varargin{1};obj=obj.PEC();obj.Type="ABZ";obj.Tensor=true;obj.Id=0;
            if(nargin==1),obj.Param=eye(3,3);
            elseif(nargin==2),param=varargin{2};
                if(isa(param,"MultiFrequency")),obj.Dispersive=true;obj.Param=cell(param.NF,1);
                    for ii=1:param.NF,obj.Param{ii}=eye(3,3);end
                elseif(iscell(param)),obj.Param=param;obj.Dispersive=true;
                else,obj.Param=param;obj.Dispersive=false;
                end
            end
        end
        function obj = Port(varargin)
            if(nargin==1),obj=varargin{1};obj=obj.PEC();obj.Type="POR";obj.Id=0;
            elseif(nargin==2),obj=varargin{1};obj.Id=0;obj.Type="POR";param=varargin{2};obj.Dispersive=true;obj.Param=zeros(param.NF,1);obj.Vec=cell(param.NF,1);
            end
        end
        function obj = Dirichlet(varargin)
            if(nargin==1),obj=varargin{1};obj=obj.PEC();obj.Type="DIR";obj.Id=1;
            elseif(nargin==2),obj=varargin{1};obj.Id=1;param=varargin{2};obj.Dispersive=true;obj.Param=zeros(param.NF,1);obj.Vec=cell(param.NF,1);obj.Type="DIR";obj.Id=1;
            end
        end
        function obj = PEC_EH(obj),obj.Dispersive=false;obj.Tensor=false;obj.Param=[];obj.Type="PEC";obj.Id=2;end
        function obj = PMC_EH(obj),obj.Dispersive=false;obj.Tensor=false;obj.Param=[];obj.Type="PMC";obj.Id=3;end
        function obj = ABC(obj),obj.Dispersive=false;obj.Tensor=false;obj.Param=[];obj.Type="ABC";obj.Id=0;end
        function obj = PBC(obj,type),obj.Dispersive=false;obj.Tensor=false;obj.Param=[];obj.Type="PBC";obj.Id=0;
            if (type=="m"),obj.Master=true;else,obj.Master=false;end
        end
        function obj = Tensorial(obj),obj.Tensor=true;
                if(obj.Dispersive),for ii=1:numel(obj.Param),obj.Param{ii}=eye(3,3);end
                else,obj.Param=eye(3,3);
                end
                if(obj.Type=="POR" || obj.Type=="DIR"),obj.PortParamType=2;end
        end
        function obj = Scalar(obj),obj.Tensor=false;
                if(obj.Dispersive),for ii=1:numel(obj.Param),obj.Param(ii)=0;end
                else,obj.Param=0;
                end
                if(obj.Type=="POR" || obj.Type=="DIR"),obj.PortParamType=1;end
        end
        function obj = Disp(obj,FRange)
                       if(obj.Tensor)
                           obj.Param=cell(FRange.NF,1);
                           for ii=1:FRange.NF,obj.Param{ii}=eye(3,3);end
                       else
                           obj.Param=zeros(FRange.NF,1);
                       end
        end
        function obj = NonDisp(obj)
                        if(obj.Tensor)
                            obj.Param=eye(3,3);
                        else
                            obj.Param=0;
                       end
        end
    end
end