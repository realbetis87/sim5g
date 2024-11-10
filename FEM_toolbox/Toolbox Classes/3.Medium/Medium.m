classdef Medium
%==========================================================================
%{
    Properties:
    
      1.Type                : "Iso"  - Isotropic Medium
                              "Anis" - Anisotropic Medium
                              "Bian" - Bi-Anisotropic Medium
      2.IsDispersive        : true/false
      3.FRange              : IsDispersive == false - null
                              IsDispersive == true  - Frequency Range vector in Hz. 
      4.Epsilon             : Relative Dielectric Permittivity er
                                  a. IsDispersive == false && Type == "Iso"
                                     scalar
                                  b. IsDispersive == false && (Type =="Anis" || Type=="Bian")
                                     (3x3) Matrix
                                  c. IsDispersive == true && Type =="Iso"
                                     (numel(FRange)x1) vector
                                  d. IsDispersive == true && (Type =="Anis" || Type=="Bian")
                                     ((numel(FRange)x1) cell
       5.Mu                 : Relative Magnetic Permeability mr
                                  a. IsDispersive == false && Type == "Iso"
                                     scalar
                                  b. IsDispersive == false && (Type =="Anis" || Type=="Bian")
                                     (3x3) Matrix
                                  c. IsDispersive == true && Type =="Iso"
                                     (numel(FRange)x1) vector
                                  d. IsDispersive == true && (Type =="Anis" || Type=="Bian")
                                     ((numel(FRange)x1) cell
      6.Ksi                 : Magneto Electric Coupling Parameter ξ
      7.Zita                : Magneto Electric Coupling Parameter ζ
      8.Tag                 : (Optional) Material Name
      9.WaveImpedance       : Medium WaveImpedance (scalar if medium.Type == "Iso")
                                                          3x3 tensor "Anis","Bian"
%}
%==========================================================================
    properties
       Type="Anis";IsDispersive=false;FRange;Epsilon=eye(3,3);Mu=eye(3,3);Ksi=zeros(3,3);Zita=zeros(3,3);Tag;WaveImpedance;
    end
    
    methods
        function obj = Medium(varargin)
            if(nargin==0),obj.Tag="Empty Medium";return;
            elseif(nargin==1),obj.Tag=varargin{1};
            elseif(nargin==2),obj.Tag=varargin{1};obj.Type=varargin{2};
            end
        end
        function obj = Dispersive(obj,FRange),obj.FRange=FRange;obj.IsDispersive=true;ElectromagneticConstants;
            if(obj.Type=="Iso"),obj.Epsilon=ones(obj.FRange.NF,1);obj.Mu=ones(obj.FRange.NF,1);obj.WaveImpedance=sqrt(m0/e0)*ones(obj.FRange.NF,1);
            elseif(obj.Type=="Anis"),obj.Epsilon=cell(obj.FRange.NF,1);obj.Mu=cell(obj.FRange.NF,1);obj.WaveImpedance=cell(obj.FRange.NF,1);
                for ii=1:obj.FRange.NF,obj.Epsilon{ii}=eye(3,3);obj.Mu{ii}=eye(3,3);obj.WaveImpedance{ii}=sqrt(m0/e0)*eye(3,3);end
            elseif(obj.Type=="Bian"),obj.Epsilon=cell(obj.FRange.NF,1);obj.Mu=cell(obj.FRange.NF,1);obj.Ksi=cell(obj.FRange.NF,1);
                obj.Zita=cell(obj.FRange.NF,1);obj.WaveImpedance=cell(obj.FRange.NF,1);
                for ii=1:obj.FRange.NF,obj.Epsilon{ii}=eye(3,3);obj.Mu{ii}=eye(3,3);obj.Ksi{ii}=zeros(3,3);obj.Zita{ii}=zeros(3,3);obj.WaveImpedance{ii}=sqrt(m0/e0)*eye(3,3);end
            end
        end
        function obj = NoDispersive(obj),obj.FRange=[];obj.IsDispersive=false;ElectromagneticConstants;
            if(obj.Type=="Iso"),obj.Epsilon=1;obj.Mu=1;obj.WaveImpedance=sqrt(m0/e0);
            elseif(obj.Type=="Anis"),obj.Epsilon=eye(3,3);obj.Mu=eye(3,3);obj.WaveImpedance=sqrt(m0/e0)*eye(3,3);
            elseif(obj.Type=="Bian"),obj.Epsilon=eye(3,3);obj.Mu=eye(3,3);obj.Ksi=zeros(3,3);obj.Zita=zeros(3,3);obj.WaveImpedance=sqrt(m0/e0)*eye(3,3);
            end
        end
        function obj = Bianisotropic(obj),obj.Type="Bian";obj.Zita=zeros(3,3);obj.Ksi=zeros(3,3);
            if(obj.IsDispersive),obj=obj.Dispersive(obj.FRange);
            else,obj=obj.NoDispersive();
            end
        end
        function obj = Isotropic(obj),obj.Type="Iso";obj.Zita=0;obj.Ksi=0;
            if(obj.IsDispersive),obj=obj.Dispersive(obj.FRange);
            else,obj=obj.NoDispersive();
            end
        end
        function obj =Anisotropic(obj),obj.Type="Anis";
            if(obj.IsDispersive),obj=obj.Dispersive(obj.FRange);
            else,obj=obj.NoDispersive();
            end
        end
        function obj =CalculateMaterialParameters(obj),ElectromagneticConstants;
            if(obj.Type=="Iso")
                if(obj.IsDispersive),for ii=1:obj.FRange.NF,obj.WaveImpedance(ii)=sqrt(obj.Mu(ii)*m0/(obj.Epsilon(ii)*e0));end
                else,obj.WaveImpedance=sqrt(obj.Mu*m0/(obj.Epsilon*e0));
                end
            elseif(obj.Type=="Anis")
                if(obj.IsDispersive),for ii=1:obj.FRange.NF,ie=e0*obj.Epsilon{ii};ie=ie^-1;sie=sqrt(ie);sm=sqrt(obj.Mu{ii}*m0);obj.WaveImpedance{ii}=sm*sie;end
                else,ie=e0*obj.Epsilon;ie=ie^-1;sie=sqrt(ie);sm=sqrt(obj.Mu*m0);obj.WaveImpedance=sm*sie;
                end
            else
            end
        end
    end
end



