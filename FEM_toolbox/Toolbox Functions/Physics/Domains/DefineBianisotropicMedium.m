%--------------------------------------------------------------------------
%{
                DefineBianisotropicMedium   
    Add isotropic medium to the toolboxModel domain with domainIndex
    
    ToolboxModel   : (ToolboxModel) object
    DomainIndex    : (int) 
    FrequencyRange : (MultiFrequency) object
    epsilon        : (double) 3x3 matrix /((FrequencyRange.NF x 1) cell)
    mu             : (double) 3x3 matrix /((FrequencyRange.NF x 1) cell)
    ksi            : (double) 3x3 matrix /((FrequencyRange.NF x 1) cell)
    zita           : (double) 3x3 matrix /((FrequencyRange.NF x 1) cell)

   1. DefineBianisotropicMedium(ToolboxModel,DomainIndex)
   2. DefineBianisotropicMedium(ToolboxModel,DomainIndex,FrequencyRange)
   3. DefineBianisotropicMedium(ToolboxModel,DomainIndex,epsilon,mu,ksi,zita)
   4. DefineBianisotropicMedium(ToolboxModel,DomainIndex,epsilon,mu,ksi,zita,FrequencyRange)
%}
%--------------------------------------------------------------------------
function [toolboxModel] = DefineBianisotropicMedium(varargin)
    if(nargin==2),toolboxModel=varargin{1};domainIndex=varargin{2};domain=toolboxModel.Domains(domainIndex);
        if(isempty(domain.Medium)),med=Medium();else,med=domain.Medium;end
        med.Bianisotropic();domain.Medium=med;toolboxModel.Domains(domainIndex)=domain;
    elseif(nargin==3),toolboxModel=varargin{1};domainIndex=varargin{2};FreqRange=varargin{3};domain=toolboxModel.Domains(domainIndex);
        if(isempty(domain.Medium)),med=Medium();else,med=domain.Medium;end
        med.Bianisotropic();med.Dispersive(FreqRange);domain.Medium=med;toolboxModel.Domains(domainIndex)=domain;
    elseif(nargin==6),toolboxModel=varargin{1};domainIndex=varargin{2};epsilon=varargin{3};mu=varargin{4};ksi=varargin{5};zita=varargin{6};
        domain=toolboxModel.Domains(domainIndex);if(isempty(domain.Medium)),med=Medium();else,med=domain.Medium;end
         med.Bianisotropic();med.Epsilon=epsilon;med.Mu=mu;med.Ksi=ksi;med.Zita=zita;domain.Medium=med;toolboxModel.Domains(domainIndex)=domain;
    elseif(nargin==7),toolboxModel=varargin{1};domainIndex=varargin{2};epsilon=varargin{3};mu=varargin{4};ksi=varargin{5};zita=varargin{6};FreqRange=varargin{7};domain=toolboxModel.Domains(domainIndex);
        if(isempty(domain.Medium)),med=Medium();else,med=domain.Medium;end
        med.Bianisotropic();med.Dispersive(FreqRange);med.Epsilon=epsilon;med.Mu=mu;med.Ksi=ksi;med.Zita=zita;       
        domain.Medium=med;toolboxModel.Domains(domainIndex)=domain;
    end
end

