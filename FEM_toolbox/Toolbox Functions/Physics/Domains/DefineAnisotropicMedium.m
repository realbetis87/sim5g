%--------------------------------------------------------------------------
%{
                DefineAnisotropicMedium   
    Add isotropic medium to the toolboxModel domain with domainIndex
    
    ToolboxModel   : (ToolboxModel) object
    DomainIndex    : (int) 
    FrequencyRange : (MultiFrequency) object
    epsilon        : (double) 3x3 matrix /((FrequencyRange.NF x 1) cell)
    mu             : (double) 3x3 matrix /((FrequencyRange.NF x 1) cell)

   1. DefineAnisotropicMedium(ToolboxModel,DomainIndex)
   2. DefineAnisotropicMedium(ToolboxModel,DomainIndex,FrequencyRange)
   3. DefineAnisotropicMedium(ToolboxModel,DomainIndex,epsilon,mu)
   4. DefineIsotropicMedium(ToolboxModel,DomainIndex,epsilon,mu,FrequencyRange)
%}
%--------------------------------------------------------------------------
function [toolboxModel] = DefineAnisotropicMedium(varargin)
    if(nargin==2),toolboxModel=varargin{1};domainIndex=varargin{2};domain=toolboxModel.Domains(domainIndex);
        if(isempty(domain.Medium)),med=Medium();else,med=domain.Medium;end
        med.Anisotropic();domain.Medium=med;toolboxModel.Domains(domainIndex)=domain;
    elseif(nargin==3),toolboxModel=varargin{1};domainIndex=varargin{2};FreqRange=varargin{3};domain=toolboxModel.Domains(domainIndex);
        if(isempty(domain.Medium)),med=Medium();else,med=domain.Medium;end
        med.Anisotropic();med.Dispersive(FreqRange);domain.Medium=med;toolboxModel.Domains(domainIndex)=domain;
    elseif(nargin==4),toolboxModel=varargin{1};domainIndex=varargin{2};epsilon=varargin{3};mu=varargin{4};domain=toolboxModel.Domains(domainIndex);
        if(isempty(domain.Medium)),med=Medium();else,med=domain.Medium;end
         med.Anisotropic();med.Epsilon=epsilon;med.Mu=mu;domain.Medium=med;toolboxModel.Domains(domainIndex)=domain;
    elseif(nargin==5),toolboxModel=varargin{1};domainIndex=varargin{2};epsilon=varargin{3};mu=varargin{4};FreqRange=varargin{5};domain=toolboxModel.Domains(domainIndex);
        if(isempty(domain.Medium)),med=Medium();else,med=domain.Medium;end
        med.Anisotropic();med.Dispersive(FreqRange);med.Epsilon=epsilon;med.Mu=mu;       
        domain.Medium=med;toolboxModel.Domains(domainIndex)=domain;
    end
end

