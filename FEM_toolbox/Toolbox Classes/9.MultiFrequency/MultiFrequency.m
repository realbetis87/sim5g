classdef MultiFrequency
%==========================================================================
%{
    1. Start            : Frequency Range Start 
    2. Stop             : Frequency Range Stop
    3. Increment        : Frequency Range Increment  
    4. Frequency        : (vector) Frequency Range (Hz)
    5. UFrequency       : (vector) Frequency Range 
    6. Unit             : string - "KHz","MHz","GHz","THz" 
    7. Omega            : (vector) angular frequency
    8. NF               : numel(Frequency)    

    Constructors 
        i.  MultiFrequency(Frequencies,unit)
        ii. MultiFrequency(start,increment,stop,unit)
%}
%==========================================================================
    properties
        Start;Stop;Increment;Frequency;Unit;Omega;NF;SingleFrequency=false;UFrequency;
    end
    
    methods
        function obj = MultiFrequency(varargin),ElectromagneticConstants;
            if (nargin==2),obj.Frequency=varargin{1};obj.Unit=varargin{2};obj.NF=numel(obj.Frequency);obj.UFrequency=varargin{1};
                if(obj.Unit=="Hz"),obj.Frequency=obj.Frequency;obj.Omega=2*pi*obj.Frequency;
                elseif(obj.Unit=="KHz"),obj.Frequency=obj.Frequency*1e3;obj.Omega=2*pi*obj.Frequency;
                elseif(obj.Unit=="MHz"),obj.Frequency=obj.Frequency*1e6;obj.Omega=2*pi*obj.Frequency;
                elseif(obj.Unit=="GHz"),obj.Frequency=obj.Frequency*1e9;obj.Omega=2*pi*obj.Frequency;
                elseif(obj.Unit=="THz"),obj.Frequency=obj.Frequency*1e12;obj.Omega=2*pi*obj.Frequency;
                end
            elseif(nargin==4),start=varargin{1};increment=varargin{2};stop=varargin{3};unit=varargin{4};obj.Start=start;obj.Stop=stop;obj.Increment=increment;obj.Frequency=start:increment:stop;obj.NF=numel(obj.Frequency);obj.Unit=unit;obj.UFrequency=obj.Frequency;
                if(unit=="Hz"),obj.Frequency=obj.Frequency;obj.Omega=2*pi*obj.Frequency;
                elseif(unit=="KHz"),obj.Frequency=obj.Frequency*1e3;obj.Omega=2*pi*obj.Frequency;
                elseif(unit=="MHz"),obj.Frequency=obj.Frequency*1e6;obj.Omega=2*pi*obj.Frequency;
                elseif(unit=="GHz"),obj.Frequency=obj.Frequency*1e9;obj.Omega=2*pi*obj.Frequency;
                elseif(unit=="THz"),obj.Frequency=obj.Frequency*1e12;obj.Omega=2*pi*obj.Frequency;
                end
            end
            if(obj.NF==1),obj.SingleFrequency=true;else,obj.SingleFrequency=false;end
        end
    end
end

