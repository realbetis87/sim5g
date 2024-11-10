classdef FEMAssembly
%==========================================================================
%{  
    Finite Element Assembled Matrices

    Properties:

    1. Type             : "Excitation"
                          "EigenFrequency"
                          "EigenMode"
    2. Dimension        : "3D","2D"
    3. Matrix_A         : Primary Matrix 
    4. Matrix_B         : Secondary Matrix
    5. NF               : NumberOfFrequencies
    6. Bian             : true/false 
    7. Dispersive       : true/false
    %----------------------- 3D Formulations ------------------------------
    1. PropagationAxis  :  "x","y","z" (EigenMode 3D)
    2. FM,AM            :  Volume Matrices (Matrix A) 
    3. TE,TB            :  Volume Matrices (Matrix A - Excitation,EigenMode)
                                           (Matrix B - EigenFrequency)   
    4. TS,TG,TBC        :  Surface Matrices (Matrix A)
    5. AW,FW            :  Volume Matrices (Matrix B) (EigenMode Formulation)
    6. P,TA,TC          :  Volume Matrices (Matrix A) (Bianisotropic Media)
    7. K                :  Volume Matrix (Matrix B) (Bianisotropic Media)
    8. TSB              :  Surface Matrix (Matrix A) (Bianisotropic Media)
    9. TP               :  Surface Matrix (Matrix A) (Port)
    10.TPV              :  Surface Matrix (Matrix B) (Port)
    11.EigenValue       : "k" - eigenvalue : Bloch-Floquet wavevector (EigenMode Formulation)
                          "n" - eigenvalue : equivalent effective
                                            refractive index
    %------------------- Dirichlet Excitation -----------------------------
    1.TEB,TBB,AB,FB     : Volume Matrices (Matrix B)  DoF Rows - Known Columns  
    2.PB,TAB,TCB        : Volume Matrices (Matrix B)  DoF Rows - Known Columns  
    %------------------- 2D EigenMode Properties --------------------------
    1.  Equation        : (2D EigenMode Formulations)                        
                            1 - 2D 1/2 EH 
                            2 - 2D EH 
                            3 - 2D 1/2 Vectorial Wave Equation for E.
    2. DimEt            : (2D FEMAssembly) number of tangential Electric DoFs
    3. DimEn            : (2D FEMAssembly) number of normal Electric DoFs
    4. DimHt            : (2D FEMAssembly) number of tangential Magnetic DoFs
    5. DimHn            : (2D FEMAssembly) number of normal Magnetic DoFs
    %----------------- Excitation Formulation -----------------------------
    1.FM,AM,TE,TB       : Volume Matrices
    2.Matrix_B          : for Dirichlet Excitation
    3.Vector            : 
    4.IsDir             : true/false
    5.IsPort            : true/false
    %----------------- Scaling Properties ---------------------------------
    9. E_Scaling       :   0 - none
                            1 - EdgeLengths
                            2 - Scalling Vector E
    10.ScalingVec_E    :   Scalling vector for E basis functions
    11.B_Scaling       :   0 - none
                            1 - FacetSurfaces
                            2 - Scalling Vector B
    12.ScalingVec_B    :   Scalling vector for B basis functions
    11. EqScalling      :  Equation Scalling (3D Assembly)
                            0  - none
                            1  - E speed of light Scalling
                            2  - B speed of light Scalling
    Functions :
    
    FEMAssembly(Type,Dimension) 
    FEMAssembly(Type,Dimension,Equation)
    FEMAssembly(Type,Dimension,FrequencyVector)
%}  
%==========================================================================
    properties
        Type;Dimension;Matrix_A;Matrix_B;NF;Bian=false;Dispersive=false;
        PropagationAxis;
        TE;TB;AM;FM;TS;TG;TBC;
        TEB;TBB;AB;FB;PB;TAB;TCB;
        TP;TPV;
        AW;FW;
        P;TA;TC;K;TSB;
        E_Scaling=1;B_Scaling=1;ScalingVec_E;ScalingVec_B;EqScaling=2;
        Equation;DimEt=0;DimHt=0;DimEn=0;DimHn=0;
        Vector;IsDir;IsPort;
        NE;NB;N;KNE;KNB;KN;
        EigenValue="k";
    end
    
    methods
        function obj = FEMAssembly(varargin)
            if(nargin==0)
            elseif(nargin==2),obj.Type=varargin{1};obj.Dimension=varargin{2};
            elseif(nargin==3),obj.Type=varargin{1};obj.Dimension=varargin{2};
                if(obj.Dimension=="3D"),obj.NF=varargin{3};
                else,obj.Equation=varargin{3};
                end
            elseif(nargin==4),obj.Type=varargin{1};obj.Dimension=varargin{2};obj.Equation=varargin{3};obj.NF=varargin{3};
            end
        end
        %---------------- Internal Functions ------------------------------
        function obj = Init3DEigenMode_MF(obj)
            obj.TE=cell(obj.NF,1);      obj.TS=cell(obj.NF,1);          obj.AW=cell(obj.NF,1);      obj.Matrix_A=cell(obj.NF,1);
            obj.TB=cell(obj.NF,1);      obj.TG=cell(obj.NF,1);          obj.FW=cell(obj.NF,1);      obj.Matrix_B=cell(obj.NF,1);
            obj.AM=cell(obj.NF,1);      obj.TBC=cell(obj.NF,1);
            obj.FM=cell(obj.NF,1);      obj.TBC=cell(obj.NF,1);
        end
        function obj = Init3DEigenModeBian_MF(obj)
            obj.TE=cell(obj.NF,1);       obj.TG=cell(obj.NF,1);         obj.P=cell(obj.NF,1);       obj.AW=cell(obj.NF,1);        obj.Matrix_A=cell(obj.NF,1);
            obj.TB=cell(obj.NF,1);       obj.TS=cell(obj.NF,1);         obj.TC=cell(obj.NF,1);      obj.FW=cell(obj.NF,1);        obj.Matrix_B=cell(obj.NF,1);
            obj.AM=cell(obj.NF,1);       obj.TBC=cell(obj.NF,1);        obj.TA=cell(obj.NF,1);      obj.K=cell(obj.NF,1);
            obj.FM=cell(obj.NF,1);       obj.TBC=cell(obj.NF,1);                                 
            obj.TSB=cell(obj.NF,1);
        end
        function obj = Init3DExcitation_MF(obj)
            if(obj.Bian==true)
                if(obj.IsDir)
                    obj.TE=cell(obj.NF,1);      obj.TS=cell(obj.NF,1);      obj.P=cell(obj.NF,1);       obj.TEB=cell(obj.NF,1);         obj.PB=cell(obj.NF,1);          obj.Matrix_A=cell(obj.NF,1);
                    obj.TB=cell(obj.NF,1);      obj.TG=cell(obj.NF,1);      obj.TA=cell(obj.NF,1);      obj.TBB=cell(obj.NF,1);         obj.TAB=cell(obj.NF,1);         obj.Matrix_B=cell(obj.NF,1);
                    obj.AM=cell(obj.NF,1);      obj.TBC=cell(obj.NF,1);     obj.TC=cell(obj.NF,1);      obj.AB=cell(obj.NF,1);          obj.TCB=cell(obj.NF,1); 
                    obj.FM=cell(obj.NF,1);      obj.TP=cell(obj.NF,1);                                  obj.FB=cell(obj.NF,1);
                                                obj.TPV=cell(obj.NF,1);
                else
                    obj.TE=cell(obj.NF,1);      obj.TS=cell(obj.NF,1);      obj.P=cell(obj.NF,1);       obj.Matrix_A=cell(obj.NF,1);
                    obj.TB=cell(obj.NF,1);      obj.TG=cell(obj.NF,1);      obj.TA=cell(obj.NF,1);      obj.Matrix_B=cell(obj.NF,1);
                    obj.AM=cell(obj.NF,1);      obj.TBC=cell(obj.NF,1);     obj.TC=cell(obj.NF,1);       
                    obj.FM=cell(obj.NF,1);      obj.TP=cell(obj.NF,1);                                  
                                                obj.TPV=cell(obj.NF,1);
                end
            else
                if(obj.IsDir)
                    obj.TE=cell(obj.NF,1);      obj.TS=cell(obj.NF,1);      obj.TEB=cell(obj.NF,1);      obj.Matrix_A=cell(obj.NF,1);
                    obj.TB=cell(obj.NF,1);      obj.TG=cell(obj.NF,1);      obj.TBB=cell(obj.NF,1);      obj.Matrix_B=cell(obj.NF,1);
                    obj.AM=cell(obj.NF,1);      obj.TBC=cell(obj.NF,1);     obj.AB=cell(obj.NF,1);          
                    obj.FM=cell(obj.NF,1);      obj.TP=cell(obj.NF,1);      obj.FB=cell(obj.NF,1);
                                                obj.TPV=cell(obj.NF,1);
                else
                    obj.TE=cell(obj.NF,1);      obj.TS=cell(obj.NF,1);      obj.Matrix_A=cell(obj.NF,1);
                    obj.TB=cell(obj.NF,1);      obj.TG=cell(obj.NF,1);      obj.Matrix_B=cell(obj.NF,1);
                    obj.AM=cell(obj.NF,1);      obj.TBC=cell(obj.NF,1);      
                    obj.FM=cell(obj.NF,1);      obj.TP=cell(obj.NF,1);                                  
                                                obj.TPV=cell(obj.NF,1);
                end

            end
        end
        function obj = Init3DEigenFrequency_MF(obj)
            obj.TE=cell(obj.NF,1);
            obj.TB=cell(obj.NF,1);obj.TG=cell(obj.NF,1);obj.TS=cell(obj.NF,1);obj.AM=cell(obj.NF,1);obj.FM=cell(obj.NF,1);
            obj.Matrix_A=cell(obj.NF,1);obj.Matrix_B=cell(obj.NF,1);
        end
        %---------------- External Functions -------------------------------
        function obj = EigenMode(varargin),ElectromagneticConstants;
                 if(nargin==12),obj=varargin{1};
                           obj.TE=varargin{2};      obj.TS=varargin{6};     obj.AW=varargin{9};
                           obj.TB=varargin{3};      obj.TG=varargin{7};     obj.FW=varargin{10};
                           obj.AM=varargin{4};      obj.TBC=varargin{8};    TModel=varargin{11};
                           obj.FM=varargin{5};                              freq=varargin{12};
                           omega=2*pi*freq;k0=omega/c0;
                           switch obj.EqScaling
                               case 0
                               case 1
                               case 2
                                   obj.Matrix_A=k0*obj.TE+1i*obj.TS+obj.FM+obj.AM+k0*obj.TB-1i*Z0*obj.TG-1i*Z0*obj.TBC;
                                   if(obj.EigenValue=="k"),obj.Matrix_B=1i*obj.FW+1i*obj.AW;
                                   else,obj.Matrix_B=1i*(k0)*obj.FW+1i*(k0)*obj.AW;
                                   end
                           end
                 elseif(nargin==13),obj=varargin{1};                    freqIndex=varargin{13};
                                    obj.TE{freqIndex}=varargin{2};      obj.TS{freqIndex}=varargin{6};        obj.AW{freqIndex}=varargin{9};
                                    obj.TB{freqIndex}=varargin{3};      obj.TG{freqIndex}=varargin{7};        obj.FW{freqIndex}=varargin{10};
                                    obj.AM{freqIndex}=varargin{4};      obj.TBC{freqIndex}=varargin{8};       TModel=varargin{11};
                                    obj.FM{freqIndex}=varargin{5};                                            freq=varargin{12};
                                    omega=2*pi*freq;k0=omega/c0;                     
                                    switch obj.EqScaling
                                        case 0
                                        case 1
                                        case 2
                                            obj.Matrix_A{freqIndex}=k0*obj.TE{freqIndex}+1i*obj.TS{freqIndex}+obj.FM{freqIndex}+obj.AM{freqIndex}+k0*obj.TB{freqIndex}-1i*Z0*obj.TG{freqIndex}-1i*Z0*obj.TBC{freqIndex};
                                            if(obj.EigenValue=="k"),obj.Matrix_B{freqIndex}=1i*obj.FW{freqIndex}+1i*obj.AW{freqIndex};
                                            else,obj.Matrix_B{freqIndex}=1i*(k0)*obj.FW{freqIndex}+1i*(k0)*obj.AW{freqIndex};
                                            end
                                    end
                 end
        end
        function obj = BianisotropicEigenMode(varargin),ElectromagneticConstants;
            if(nargin==16),obj=varargin{1}; 
                           obj.TE=varargin{2};      obj.TS=varargin{6};     obj.AW=varargin{9};     obj.P=varargin{11};        TModel=varargin{15};
                           obj.TB=varargin{3};      obj.TG=varargin{7};     obj.FW=varargin{10};    obj.TA=varargin{12};       freq=varargin{16};
                           obj.AM=varargin{4};      obj.TBC=varargin{8};                            obj.TC=varargin{13};   
                           obj.FM=varargin{5};                                                      obj.K=varargin{14};
                           omega=2*pi*freq;k0=omega/c0;                   
                           switch obj.EqScaling
                               case 0
                               case 1
                               case 2
                                   obj.Matrix_A=k0*obj.TE+1i*obj.TS+obj.FM+obj.AM+k0*obj.TB-omega*obj.TC-omega*c0*obj.TA-1i*c0*obj.P-1i*Z0*obj.TG-1i*Z0*obj.TBC;
                                   if(obj.EigenValue=="k"),obj.Matrix_B=1i*obj.FW+1i*obj.AW+c0*obj.K;
                                   else,obj.Matrix_B=(1i*k0)*obj.FW+(1i*k0)*obj.AW+(c0*k0)*obj.K;
                                   end
                           end
            elseif(nargin==17),obj=varargin{1};                    freqIndex=varargin{17};
                               obj.TE{freqIndex}=varargin{2};      obj.TS{freqIndex}=varargin{6};     obj.AW{freqIndex}=varargin{9};     obj.P{freqIndex}=varargin{11};        TModel=varargin{15};
                               obj.TB{freqIndex}=varargin{3};      obj.TG{freqIndex}=varargin{7};     obj.FW{freqIndex}=varargin{10};    obj.TA{freqIndex}=varargin{12};       freq=varargin{16};
                               obj.AM{freqIndex}=varargin{4};      obj.TBC{freqIndex}=varargin{8};                                       obj.TC{freqIndex}=varargin{13};   
                               obj.FM{freqIndex}=varargin{5};                                                                            obj.K{freqIndex}=varargin{14};
                               omega=2*pi*freq;k0=omega/c0;
                               switch obj.EqScaling
                                        case 0
                                        case 1
                                        case 2
                                            obj.Matrix_A{freqIndex}=k0*obj.TE{freqIndex}+1i*obj.TS{freqIndex}+obj.FM{freqIndex}+obj.AM{freqIndex}+k0*obj.TB{freqIndex}-1i*Z0*obj.TG{freqIndex}-1i*Z0*obj.TBC{freqIndex};
                                            if(obj.EigenValue=="k"),obj.Matrix_B{freqIndex}=1i*obj.FW{freqIndex}+1i*obj.AW{freqIndex}+c0*obj.K{freqIndex};
                                            else,obj.Matrix_B=1i*(k0)*obj.FW{freqIndex}+1i*(k0)*obj.AW{freqIndex}+(c0*k0)*obj.K{freqIndex};
                                            end
                               end
            end
        end
        function obj = SetScalling(varargin)
            if(nargin==3),obj=varargin{1};type=varargin{2};val=varargin{3};
                switch type
                    case "Equation", obj.EqScalling=val; 
                    case "BasisFunction",if(isvector(val)),obj.BFScaling=5;obj.ScalingVec=val;else,obj.BFScaling=val;end
                end
            elseif(nargin==5),obj=varargin{1};type=varargin{2};val=varargin{3};type2=varargin{4};val2=varargin{5};
                switch type
                    case "Equation",obj.EqScaling=val;
                    case "BasisFunction",if(isvector(val)),obj.BFScaling=5;obj.ScalingVec=val;else,obj.BFScaling=val;end
                end
                switch type2
                    case "Equation",obj.EqScaling=val2;
                    case "BasisFunction",if(isvector(val2)),obj.BFScaling=5;obj.ScalingVec=val2;else,obj.BFScaling=val2;end
                end
            end
        end
        function obj = DirichletExcitation(varargin),ElectromagneticConstants;
            if(nargin==15),obj=varargin{1};
                           obj.TE=varargin{2};      obj.TS=varargin{6};     obj.TEB=varargin{10};        TModel=varargin{14};
                           obj.TB=varargin{3};      obj.TP=varargin{7};     obj.TBB=varargin{11};        frequency=varargin{15};
                           obj.AM=varargin{4};      obj.TG=varargin{8};     obj.AB=varargin{12};
                           obj.FM=varargin{5};      obj.TBC=varargin{9};    obj.FB=varargin{13};
                           
                           k0=2*pi*frequency/c0;
                           switch obj.EqScaling
                              case 0
                              case 1
                              case 2
                                  obj.Matrix_A=k0*obj.TE+1i*obj.TS+obj.FM+obj.AM+k0*obj.TB-1i*Z0*obj.TG-1i*Z0*obj.TBC;
                                  obj.TPV=obj.TP;obj.Matrix_B=k0*obj.TEB+obj.FB+obj.AB+k0*obj.TBB;
                          end
            elseif(nargin==16),obj=varargin{1};                 freqIndex=varargin{16};
                           obj.TE{freqIndex}=varargin{2};       obj.TS{freqIndex}=varargin{6};      obj.TEB{freqIndex}=varargin{10};    TModel=varargin{14};
                           obj.TB{freqIndex}=varargin{3};       obj.TP{freqIndex}=varargin{7};      obj.TBB{freqIndex}=varargin{11};    frequency=varargin{15};    
                           obj.FM{freqIndex}=varargin{4};       obj.TG{freqIndex}=varargin{8};      obj.AB{freqIndex}=varargin{12};
                           obj.AM{freqIndex}=varargin{5};       obj.TBC{freqIndex}=varargin{9};     obj.FB{freqIndex}=varargin{13};
                                         
                           k0=2*pi*frequency/c0;
                           switch obj.EqScaling
                              case 0
                              case 1
                              case 2
                                  obj.Matrix_A{freqIndex}=k0*obj.TE{freqIndex}+k0*obj.TB{freqIndex}+obj.AM{freqIndex}+obj.FM{freqIndex}-1i*Z0*obj.TBC{freqIndex}-1i*Z0*obj.TG{freqIndex}+1i*obj.TS{freqIndex};
                                  obj.TPV{freqIndex}=obj.TP{freqIndex};obj.Matrix_B{freqIndex}=k0*obj.TEB{freqIndex}+obj.FB{freqIndex}+obj.AB{freqIndex}+k0*obj.TBB{freqIndex};
                          end
            end 
            
        end
        function obj = BianisotropicDirichletExcitation(varargin),ElectromagneticConstants;
               if(nargin==21),obj=varargin{1}; 
                           obj.TE=varargin{2};      obj.TS=varargin{6};        obj.P=varargin{10};      obj.TEB=varargin{13};     obj.PB=varargin{17};      TModel=varargin{20};
                           obj.TB=varargin{3};      obj.TP=varargin{7};        obj.TA=varargin{11};     obj.TBB=varargin{14};     obj.TAB=varargin{18};     freq=varargin{21};
                           obj.AM=varargin{4};      obj.TG=varargin{8};        obj.TC=varargin{12};     obj.AB=varargin{15};      obj.TCB=varargin{19};  
                           obj.FM=varargin{5};      obj.TBC=varargin{9};                                obj.FB=varargin{16};               
                           omega=2*pi*freq;k0=omega/c0;    
                           switch obj.EqScaling
                               case 0
                               case 1
                               case 2
                                   obj.Matrix_A=k0*obj.TE+1i*obj.TS+obj.FM+obj.AM+k0*obj.TB-omega*obj.TC-omega*c0*obj.TA-1i*c0*obj.P-1i*Z0*obj.TG-1i*Z0*obj.TBC;
                                   obj.TPV=obj.TP;obj.Matrix_B=k0*obj.TEB+obj.FB+obj.AB+k0*obj.TBB-omega*obj.TCB-omega*c0*obj.TAB-1i*c0*obj.PB;
                           end
               elseif(nargin==22),obj=varargin{1};                    freqIndex=varargin{22};
                                  obj.TE{freqIndex}=varargin{2};      obj.TS{freqIndex}=varargin{6};        obj.P{freqIndex}=varargin{10};      obj.TEB{freqIndex}=varargin{13};     obj.PB{freqIndex}=varargin{17};      TModel=varargin{15};
                                  obj.TB{freqIndex}=varargin{3};      obj.TP{freqIndex}=varargin{7};        obj.TA{freqIndex}=varargin{11};     obj.TBB{freqIndex}=varargin{14};     obj.TAB{freqIndex}=varargin{18};     freq=varargin{16};
                                  obj.AM{freqIndex}=varargin{4};      obj.TG{freqIndex}=varargin{8};        obj.TC{freqIndex}=varargin{12};     obj.AB{freqIndex}=varargin{15};      obj.TCB{freqIndex}=varargin{19};  
                                  obj.FM{freqIndex}=varargin{5};      obj.TBC{freqIndex}=varargin{9};                                           obj.FB{freqIndex}=varargin{16};               
                                  omega=2*pi*freq;k0=omega/c0;    
                                  switch obj.EqScaling
                                        case 0
                                        case 1
                                        case 2
                                            obj.Matrix_A{freqIndex}=k0*obj.TE{freqIndex}+1i*obj.TS{freqIndex}+obj.FM{freqIndex}+obj.AM{freqIndex}+k0*obj.TB{freqIndex}-omega*obj.TC{freqIndex}-omega*c0*obj.TA{freqIndex}-1i*c0*obj.P{freqIndex}-1i*Z0*obj.TG{freqIndex}-1i*Z0*obj.TBC{freqIndex};
                                            obj.TPV{freqIndex}=obj.TP{freqIndex};obj.Matrix_B{freqIndex}=k0*obj.TEB{freqIndex}+obj.FB{freqIndex}+obj.AB{freqIndex}+k0*obj.TBB{freqIndex}-omega*obj.TCB{freqIndex}-omega*c0*obj.TAB{freqIndex}-1i*c0*obj.PB{freqIndex};
                                   end
               end
        end
        function obj = PortExcitation(varargin),obj=varargin{1};ElectromagneticConstants;
                    if(nargin==11),obj=varargin{1};
                                   obj.TE=varargin{2};      obj.TS=varargin{6};     TModel=varargin{10};
                                   obj.TB=varargin{3};      obj.TP=varargin{7};     frequency=varargin{11};
                                   obj.AM=varargin{4};      obj.TG=varargin{8};     
                                   obj.FM=varargin{5};      obj.TBC=varargin{9};    
                           
                           k0=2*pi*frequency/c0;
                           switch obj.EqScaling
                              case 0
                              case 1
                              case 2
                                  obj.Matrix_A=k0*obj.TE+1i*obj.TS+obj.FM+obj.AM+k0*obj.TB-1i*Z0*obj.TG-1i*Z0*obj.TBC;
                                  obj.TPV=obj.TP;
                          end
                    elseif(nargin==12),obj=varargin{1};                     freqIndex=varargin{12};
                                       obj.TE{freqIndex}=varargin{2};       obj.TS{freqIndex}=varargin{6};      TModel=varargin{10};
                                       obj.TB{freqIndex}=varargin{3};       obj.TP{freqIndex}=varargin{7};      frequency=varargin{11};
                                       obj.FM{freqIndex}=varargin{4};       obj.TG{freqIndex}=varargin{8};     
                                       obj.AM{freqIndex}=varargin{5};       obj.TBC{freqIndex}=varargin{9};
                                         
                           k0=2*pi*frequency/c0;
                           switch obj.EqScaling
                              case 0
                              case 1
                              case 2
                                  obj.Matrix_A{freqIndex}=k0*obj.TE{freqIndex}+k0*obj.TB{freqIndex}+obj.AM{freqIndex}+obj.FM{freqIndex}-1i*Z0*obj.TBC{freqIndex}-1i*Z0*obj.TG{freqIndex}+1i*obj.TS{freqIndex};
                                  obj.TPV{freqIndex}=obj.TP{freqIndex};
                          end
                    end 
        end
        function obj = BianisotropicPortExcitation(varargin)
                 if(nargin==14),obj=varargin{1};
                                obj.TE=varargin{2};      obj.TS=varargin{6};    obj.P=varargin{10};   TModel=varargin{13};
                                obj.TB=varargin{3};      obj.TP=varargin{7};    obj.TA=varargin{11};  frequency=varargin{14};
                                obj.AM=varargin{4};      obj.TG=varargin{8};    obj.TC=varargin{12}; 
                                obj.FM=varargin{5};      obj.TBC=varargin{9};    
                           
                           k0=2*pi*frequency/c0;
                           switch obj.EqScaling
                              case 0
                              case 1
                              case 2
                                  obj.Matrix_A=k0*obj.TE+1i*obj.TS+obj.FM+obj.AM+k0*obj.TB-omega*obj.TC-omega*c0*obj.TA-1i*c0*obj.P-1i*Z0*obj.TG-1i*Z0*obj.TBC;
                                  obj.TPV=obj.TP;
                          end

                 elseif(narginn==15),obj=varargin{1};         freqIndex=varargin{15};
                                     obj.TE{freqIndex}=varargin{2};      obj.TS{freqIndex}=varargin{6};    obj.P{freqIndex}=varargin{10};   TModel=varargin{13};
                                     obj.TB{freqIndex}=varargin{3};      obj.TP{freqIndex}=varargin{7};    obj.TA{freqIndex}=varargin{11};  frequency=varargin{14};
                                     obj.AM{freqIndex}=varargin{4};      obj.TG{freqIndex}=varargin{8};    obj.TC{freqIndex}=varargin{12}; 
                                     obj.FM{freqIndex}=varargin{5};      obj.TBC{freqIndex}=varargin{9};    
                           
                           k0=2*pi*frequency/c0;
                           switch obj.EqScaling
                              case 0
                              case 1
                              case 2
                                  obj.Matrix_A{freqIndex}=k0*obj.TE{freqIndex}+1i*obj.TS{freqIndex}+obj.FM{freqIndex}+obj.AM{freqIndex}+k0*obj.TB{freqIndex}-omega*obj.TC{freqIndex}-omega*c0*obj.TA{freqIndex}-1i*c0*obj.P{freqIndex}-1i*Z0*obj.TG{freqIndex}-1i*Z0*obj.TBC{freqIndex};
                                  obj.TPV{freqIndex}=obj.TP{freqIndex};
                          end
                 end
        end
    end
end

