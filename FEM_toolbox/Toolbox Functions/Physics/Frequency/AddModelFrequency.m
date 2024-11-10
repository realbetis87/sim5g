%--------------------------------------------------------------------------
%{
%}
%--------------------------------------------------------------------------
function ToolboxModel = AddModelFrequency(varargin)
    if(nargin==3),ToolboxModel=varargin{1};frequency=varargin{2};Unit=varargin{3};Frequency=MultiFrequency(frequency,Unit);ToolboxModel=ToolboxModel.AddFrequency(Frequency);
    elseif(nargin==5),ToolboxModel=varargin{1};Start=varargin{2};Increment=varargin{3};Stop=varargin{4};Unit=varargin{5};Frequency=MultiFrequency(Start,Increment,Stop,Unit);ToolboxModel=ToolboxModel.AddFrequency(Frequency);
    end
end