classdef Solution
%==========================================================================
%{
    1. Type           :  "EigenValue","EigenMode","Excitation"
    2. FRange         :   MultiFrequency Object (Redundant) 
    %------------------- EigenMode-----------------------------------------
    1.EigenVectors    : Single Frequency Operation : (Nv x Nu) matrix 
                        Nv number of EigenVectors , Nu Number of Degrees of Fredom
                        Frequency Range Operation : Nf x1 cell 
                        (Nf Number of Frequencies in frequency range)
    2.EigenValues
    3.Type            : "k" (eigevalue is wavevector)
                        "n" (eigenvalue is effective refractive index)
    4.Shifts
    %------------------- Excitation----------------------------------------
    1.KnownExcition
    2.UknownExcitation
    2.SolutionVector
%}
%==========================================================================
    properties
        Type;EigenValues;EigenVectors;FRange;KnownExcitation;UknownExcitation;SolutionVector;ExcitationVector;EigenValueType;
    end
    
    methods
        function obj = Solution()
        end
    end
end

