%--------------------------------------------------------------------------
%{
                DefineIsotropicMedium   
    Add isotropic medium to the toolboxModel domain with domainIndex
    
    ToolboxModel   : (ToolboxModel) object
    DomainIndex    : (int) 
    FrequencyRange : (MultiFrequency) object
    epsilon        : (double)/((FrequencyRange.NF x 1) vector)
    mu             : (double)/((FrequencyRange.NF x 1) vector)

   1. DefineIsotropicMedium(ToolboxModel,DomainIndex)
   2. DefineIsotropicMedium(ToolboxModel,DomainIndex,FrequencyRange)
   3. DefineIsotropicMedium(ToolboxModel,DomainIndex,epsilon,mu)
   4. DefineIsotropicMedium(ToolboxModel,DomainIndex,epsilon,mu,FrequencyRange)
%}
%--------------------------------------------------------------------------
function [toolboxModel] = DefineIsotropicMedium(varargin)
    if(nargin==2),toolboxModel=varargin{1};domainIndex=varargin{2};domain=toolboxModel.Domains(domainIndex);
        if(isempty(domain.Medium)),med=Medium();else,med=domain.Medium;end
        med.Isotropic();domain.Medium=med;toolboxModel.Domains(domainIndex)=domain;
    elseif(nargin==3),toolboxModel=varargin{1};domainIndex=varargin{2};FreqRange=varargin{3};domain=toolboxModel.Domains(domainIndex);
        if(isempty(domain.Medium)),med=Medium();else,med=domain.Medium;end
        med.Isotropic();med.Dispersive(FreqRange);domain.Medium=med;toolboxModel.Domains(domainIndex)=domain;
    elseif(nargin==4),toolboxModel=varargin{1};domainIndex=varargin{2};epsilon=varargin{3};mu=varargin{4};domain=toolboxModel.Domains(domainIndex);
        if(isempty(domain.Medium)),med=Medium();else,med=domain.Medium;end
         med.Isotropic();med.Epsilon=epsilon;med.Mu=mu;domain.Medium=med;toolboxModel.Domains(domainIndex)=domain;
    elseif(nargin==5),toolboxModel=varargin{1};domainIndex=varargin{2};epsilon=varargin{3};mu=varargin{4};FreqRange=varargin{5};domain=toolboxModel.Domains(domainIndex);
        if(isempty(domain.Medium)),med=Medium();else,med=domain.Medium;end
        med.Isotropic();med.Dispersive(FreqRange);med.Epsilon=epsilon;med.Mu=mu;       
        domain.Medium=med;toolboxModel.Domains(domainIndex)=domain;
    end
end

