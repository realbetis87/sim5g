%--------------------------------------------------------------------------
%{
            2D 1/2 Vector Wave Equation in terms of E
    
       Assembly for the 2D 1/2 Vector Wave Equation in Terms of the intensity 
       of the electric field E. Corresponding Problem is a linearized
       Quadratic Eigenvalue Problem :
    
        ((lamda^2)*[M] +lamda*[L]+[K])*X = 0 (X=[Et;En];

       -> Bianisotropic Media : BVWE_Assembly_Boundary
       -> Anisotropic   Media : AVWE_Assembly_Boundary
       -> Isotropic     Media : IVWE_Assembly_Boundary   
%}
%--------------------------------------------------------------------------
function [boundaryExcitation] = VWE2D_Assembly(toolboxModel,boundaryExcitation,boundaryIndices)
    if(boundaryExcitation.Dispersive == false)
        boundary=toolboxModel.Boundaries(boundaryIndices(1));
        switch abs(boundary.Axis)
           case 1,boundaryExcitation.Assembly = VWE_Assembly_xBoundary(toolboxModel,boundaryExcitation.Assembly,boundaryIndices,sign(boundary.Axis));
           case 2,boundaryExcitation.Assembly = VWE_Assembly_yBoundary(toolboxModel,boundaryExcitation.Assembly,boundaryIndices,sign(boundary.Axis));
           case 3,boundaryExcitation.Assembly = VWE_Assembly_zBoundary(toolboxModel,boundaryExcitation.Assembly,boundaryIndices,sign(boundary.Axis));
        end
    else
        boundary=toolboxModel.Boundaries(boundaryIndices(1));
        switch abs(boundary.Axis)
            case 1,for ii=1:toolboxModel.Frequency.NF,boundaryExcitation.Assembly = VWE_Assembly_xBoundary(toolboxModel,boundaryExcitation.Assembly,boundaryIndices,sign(boundary.Axis),ii);end
            case 2,for ii=1:toolboxModel.Frequency.NF,boundaryExcitation.Assembly = VWE_Assembly_yBoundary(toolboxModel,boundaryExcitation.Assembly,boundaryIndices,sign(boundary.Axis),ii);end
            case 3,for ii=1:toolboxModel.Frequency.NF,boundaryExcitation.Assembly = VWE_Assembly_zBoundary(toolboxModel,boundaryExcitation.Assembly,boundaryIndices,sign(boundary.Axis),ii);end
        end
    end
end
%===================== Boundary Assemblies ================================
function Assembly = VWE_Assembly_xBoundary(varargin)
    if(nargin==4)
        toolboxModel=varargin{1};Assembly=varargin{2};boundaryIndices=varargin{3};axisSign=varargin{4};N=Assembly.DimEt+Assembly.DimEn;

        IK=zeros(120*9*N,1);    IL=zeros(120*9*N,1);        IM=zeros(120*9*N,1);
        JK=zeros(120*9*N,1);    JL=zeros(120*9*N,1);        JM=zeros(120*9*N,1);
        K=zeros(120*9*N,1);     L=zeros(120*9*N,1);         M=zeros(120*9*N,1);
        counterK=0;             counterL=0;                 counterM=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case "Iso" ,[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_xBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign);
                    case "Anis",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_xBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign);
                    case "Bian",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_xBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign);
                end
            end
        end
        Nv=Assembly.DimEn;
        IK(counterK+1)=Nv;              IL(counterL+1)=Nv;              IM(counterM+1)=Nv;
        JK(counterK+1)=Nv;              JL(counterL+1)=Nv;              JM(counterM+1)=Nv;
        K(counterK+1)=0;                L(counterL+1)=Nv;               M(counterM+1)=0;
            
        nonZerosK=nnz(IK);              nonZerosL=nnz(IL);               nonZerosM=nnz(IM);
        IK=IK(1:nonZerosK);             IL=IL(1:nonZerosL);              IM=IM(1:nonZerosM);
        JK=JK(1:nonZerosK);             JL=JL(1:nonZerosL);              JM=JM(1:nonZerosM);
        K=K(1:nonZerosK);               L=L(1:nonZerosL);                M=M(1:nonZerosM);
        K=sparse(IK,JK,K);              L=sparse(IL,JL,L);               M=sparse(IM,JM,M);

        SA=sparse(2*Nv,2*Nv);           SB=sparse(2*Nv,2*Nv);
        SA(1:Nv,Nv+1:end)=-K;           SB(1:Nv,1:Nv)=-K;
        SA(1+Nv:end,1:Nv)=-K;           SB(Nv+1:end,Nv+1:end)=M;
        SA(1+Nv:end,1+Nv:end)=-L;

        Assembly.Matrix_A=SA;    Assembly.Matrix_B=SB;
    elseif(nargin==5)
        toolboxModel=varargin{1};Assembly=varargin{2};boundaryIndices=varargin{3};axisSign=varargin{4};frequencyIndex=varargin{5};N=Assembly.DimEt+Assembly.DimEn;
        IK=zeros(120*9*N,1);    IL=zeros(120*9*N,1);        IM=zeros(120*9*N,1);
        JK=zeros(120*9*N,1);    JL=zeros(120*9*N,1);        JM=zeros(120*9*N,1);
        K=zeros(120*9*N,1);     L=zeros(120*9*N,1);         M=zeros(120*9*N,1);
        counterK=0;             counterL=0;                 counterM=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case "Iso" ,[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_xBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign,frequencyIndex);
                    case "Anis",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_xBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign,frequencyIndex);
                    case "Bian",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_xBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign,frequencyIndex);
                end
            end
        end
        Nv=Assembly.DimEn;
        IK(counterK+1)=Nv;              IL(counterL+1)=Nv;              IM(counterM+1)=Nv;
        JK(counterK+1)=Nv;              JL(counterL+1)=Nv;              JM(counterM+1)=Nv;
        K(counterK+1)=0;                L(counterL+1)=Nv;               M(counterM+1)=0;
            
        nonZerosK=nnz(IK);              nonZerosL=nnz(IL);               nonZerosM=nnz(IM);
        IK=IK(1:nonZerosK);             IL=IL(1:nonZerosL);              IM=IM(1:nonZerosM);
        JK=JK(1:nonZerosK);             JL=JL(1:nonZerosL);              JM=JM(1:nonZerosM);
        K=K(1:nonZerosK);               L=L(1:nonZerosL);                M=M(1:nonZerosM);
        K=sparse(IK,JK,K);              L=sparse(IL,JL,L);               M=sparse(IM,JM,M);

        SA=sparse(2*Nv,2*Nv);           SB=sparse(2*Nv,2*Nv);
        SA(1:Nv,Nv+1:end)=-K;           SB(1:Nv,1:Nv)=-K;
        SA(1+Nv:end,1:Nv)=-K;           SB(Nv+1:end,Nv+1:end)=M;
        SA(1+Nv:end,1+Nv:end)=-L;

        Assembly.Matrix_A{frequencyIndex}=SA;    Assembly.Matrix_B{frequencyIndex}=SB;
    end
end
function Assembly = VWE_Assembly_yBoundary(varargin)
    if(nargin==4)
        toolboxModel=varargin{1};Assembly=varargin{2};boundaryIndices=varargin{3};axisSign=varargin{4};N=Assembly.DimEt+Assembly.DimEn;

        IK=zeros(120*9*N,1);    IL=zeros(120*9*N,1);        IM=zeros(120*9*N,1);
        JK=zeros(120*9*N,1);    JL=zeros(120*9*N,1);        JM=zeros(120*9*N,1);
        K=zeros(120*9*N,1);     L=zeros(120*9*N,1);         M=zeros(120*9*N,1);
        counterK=0;             counterL=0;                 counterM=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case "Iso" ,[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_yBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign);
                    case "Anis",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_yBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign);
                    case "Bian",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_yBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign);
                end
            end
        end
        Nv=Assembly.DimEn;
        IK(counterK+1)=Nv;              IL(counterL+1)=Nv;              IM(counterM+1)=Nv;
        JK(counterK+1)=Nv;              JL(counterL+1)=Nv;              JM(counterM+1)=Nv;
        K(counterK+1)=0;                L(counterL+1)=Nv;               M(counterM+1)=0;
            
        nonZerosK=nnz(IK);              nonZerosL=nnz(IL);               nonZerosM=nnz(IM);
        IK=IK(1:nonZerosK);             IL=IL(1:nonZerosL);              IM=IM(1:nonZerosM);
        JK=JK(1:nonZerosK);             JL=JL(1:nonZerosL);              JM=JM(1:nonZerosM);
        K=K(1:nonZerosK);               L=L(1:nonZerosL);                M=M(1:nonZerosM);
        K=sparse(IK,JK,K);              L=sparse(IL,JL,L);               M=sparse(IM,JM,M);

        SA=sparse(2*Nv,2*Nv);           SB=sparse(2*Nv,2*Nv);
        SA(1:Nv,Nv+1:end)=-K;           SB(1:Nv,1:Nv)=-K;
        SA(1+Nv:end,1:Nv)=-K;           SB(Nv+1:end,Nv+1:end)=M;
        SA(1+Nv:end,1+Nv:end)=-L;

        Assembly.Matrix_A=SA;    Assembly.Matrix_B=SB;
    elseif(nargin==5)
        toolboxModel=varargin{1};Assembly=varargin{2};boundaryIndices=varargin{3};axisSign=varargin{4};frequencyIndex=varargin{5};N=Assembly.DimEt+Assembly.DimEn;
        IK=zeros(120*9*N,1);    IL=zeros(120*9*N,1);        IM=zeros(120*9*N,1);
        JK=zeros(120*9*N,1);    JL=zeros(120*9*N,1);        JM=zeros(120*9*N,1);
        K=zeros(120*9*N,1);     L=zeros(120*9*N,1);         M=zeros(120*9*N,1);
        counterK=0;             counterL=0;                 counterM=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case "Iso" ,[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_yBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign,frequencyIndex);
                    case "Anis",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_yBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign,frequencyIndex);
                    case "Bian",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_yBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign,frequencyIndex);
                end
            end
        end
        Nv=Assembly.DimEn;
        IK(counterK+1)=Nv;              IL(counterL+1)=Nv;              IM(counterM+1)=Nv;
        JK(counterK+1)=Nv;              JL(counterL+1)=Nv;              JM(counterM+1)=Nv;
        K(counterK+1)=0;                L(counterL+1)=Nv;               M(counterM+1)=0;
            
        nonZerosK=nnz(IK);              nonZerosL=nnz(IL);               nonZerosM=nnz(IM);
        IK=IK(1:nonZerosK);             IL=IL(1:nonZerosL);              IM=IM(1:nonZerosM);
        JK=JK(1:nonZerosK);             JL=JL(1:nonZerosL);              JM=JM(1:nonZerosM);
        K=K(1:nonZerosK);               L=L(1:nonZerosL);                M=M(1:nonZerosM);
        K=sparse(IK,JK,K);              L=sparse(IL,JL,L);               M=sparse(IM,JM,M);

        SA=sparse(2*Nv,2*Nv);           SB=sparse(2*Nv,2*Nv);
        SA(1:Nv,Nv+1:end)=-K;           SB(1:Nv,1:Nv)=-K;
        SA(1+Nv:end,1:Nv)=-K;           SB(Nv+1:end,Nv+1:end)=M;
        SA(1+Nv:end,1+Nv:end)=-L;

        Assembly.Matrix_A{frequencyIndex}=SA;    Assembly.Matrix_B{frequencyIndex}=SB;
    end
end
function Assembly = VWE_Assembly_zBoundary(varargin)
    if(nargin==4)
        toolboxModel=varargin{1};Assembly=varargin{2};boundaryIndices=varargin{3};axisSign=varargin{4};N=Assembly.DimEt+Assembly.DimEn;

        IK=zeros(120*9*N,1);    IL=zeros(120*9*N,1);        IM=zeros(120*9*N,1);
        JK=zeros(120*9*N,1);    JL=zeros(120*9*N,1);        JM=zeros(120*9*N,1);
        K=zeros(120*9*N,1);     L=zeros(120*9*N,1);         M=zeros(120*9*N,1);
        counterK=0;             counterL=0;                 counterM=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case "Iso" ,[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_zBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign);
                    case "Anis",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_zBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign);
                    case "Bian",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_zBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign);
                end
            end
        end
        Nv=Assembly.DimEn;
        IK(counterK+1)=Nv;              IL(counterL+1)=Nv;              IM(counterM+1)=Nv;
        JK(counterK+1)=Nv;              JL(counterL+1)=Nv;              JM(counterM+1)=Nv;
        K(counterK+1)=0;                L(counterL+1)=Nv;               M(counterM+1)=0;
            
        nonZerosK=nnz(IK);              nonZerosL=nnz(IL);               nonZerosM=nnz(IM);
        IK=IK(1:nonZerosK);             IL=IL(1:nonZerosL);              IM=IM(1:nonZerosM);
        JK=JK(1:nonZerosK);             JL=JL(1:nonZerosL);              JM=JM(1:nonZerosM);
        K=K(1:nonZerosK);               L=L(1:nonZerosL);                M=M(1:nonZerosM);
        K=sparse(IK,JK,K);              L=sparse(IL,JL,L);               M=sparse(IM,JM,M);

        SA=sparse(2*Nv,2*Nv);           SB=sparse(2*Nv,2*Nv);
        SA(1:Nv,Nv+1:end)=-K;           SB(1:Nv,1:Nv)=-K;
        SA(1+Nv:end,1:Nv)=-K;           SB(Nv+1:end,Nv+1:end)=M;
        SA(1+Nv:end,1+Nv:end)=-L;

        Assembly.Matrix_A=SA;    Assembly.Matrix_B=SB;
    elseif(nargin==5)
        toolboxModel=varargin{1};Assembly=varargin{2};boundaryIndices=varargin{3};axisSign=varargin{4};frequencyIndex=varargin{5};N=Assembly.DimEt+Assembly.DimEn;
        IK=zeros(120*9*N,1);    IL=zeros(120*9*N,1);        IM=zeros(120*9*N,1);
        JK=zeros(120*9*N,1);    JL=zeros(120*9*N,1);        JM=zeros(120*9*N,1);
        K=zeros(120*9*N,1);     L=zeros(120*9*N,1);         M=zeros(120*9*N,1);
        counterK=0;             counterL=0;                 counterM=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case "Iso" ,[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_zBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign,frequencyIndex);
                    case "Anis",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_zBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign,frequencyIndex);
                    case "Bian",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_zBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,axisSign,frequencyIndex);
                end
            end
        end
        Nv=Assembly.DimEn;
        IK(counterK+1)=Nv;              IL(counterL+1)=Nv;              IM(counterM+1)=Nv;
        JK(counterK+1)=Nv;              JL(counterL+1)=Nv;              JM(counterM+1)=Nv;
        K(counterK+1)=0;                L(counterL+1)=Nv;               M(counterM+1)=0;
            
        nonZerosK=nnz(IK);              nonZerosL=nnz(IL);               nonZerosM=nnz(IM);
        IK=IK(1:nonZerosK);             IL=IL(1:nonZerosL);              IM=IM(1:nonZerosM);
        JK=JK(1:nonZerosK);             JL=JL(1:nonZerosL);              JM=JM(1:nonZerosM);
        K=K(1:nonZerosK);               L=L(1:nonZerosL);                M=M(1:nonZerosM);
        K=sparse(IK,JK,K);              L=sparse(IL,JL,L);               M=sparse(IM,JM,M);

        SA=sparse(2*Nv,2*Nv);           SB=sparse(2*Nv,2*Nv);
        SA(1:Nv,Nv+1:end)=-K;           SB(1:Nv,1:Nv)=-K;
        SA(1+Nv:end,1:Nv)=-K;           SB(Nv+1:end,Nv+1:end)=M;
        SA(1+Nv:end,1+Nv:end)=-L;

        Assembly.Matrix_A{frequencyIndex}=SA;    Assembly.Matrix_B{frequencyIndex}=SB;
    end
end
%===================== Isotropic Assemblies ===============================
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM] = IVWE_Assembly_xBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==15)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 

                medium = element.Medium2D; epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;
    elseif(nargin==16)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 
                                                                                                                        freqIndex= varargin{16};
                medium = element.Medium2D; 
                if(medium.IsDispersive),epsilon = medium.Epsilon(freqIndex); mu = medium.Mu(freqIndex); waveImpedance = medium.WaveImpedance(freqIndex);
                else,epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;
                end
    end
    imu=(m0^-1)*mu^-1;epsilon=e0*epsilon;omega=2*pi*frequency;
    vertices =[toolboxModel.Vertices(element.Vertices)];
    edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];z=[vertices.Z];De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;    c(1)=(y(3)-y(2))/De;
    b(2)=(z(3)-z(1))/De;    c(2)=(y(1)-y(3))/De;
    b(3)=(z(1)-z(2))/De;    c(3)=(y(2)-y(1))/De;  
    %----------------------------------------------------------------------
    Stt=zeros(3,3);         Ttt=zeros(3,3);         Ltz=zeros(3,3);         Ktt=zeros(3,3);
    Szz=zeros(3,3);         Tzz=zeros(3,3);         Lzt=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwy(1)=simp(1)*b(2)-simp(2)*b(1);       wwz(1)=simp(1)*c(2)-simp(2)*c(1);
        wwy(2)=simp(2)*b(3)-simp(3)*b(2);       wwz(2)=simp(2)*c(3)-simp(3)*c(2);
        wwy(3)=simp(3)*b(1)-simp(1)*b(3);       wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);      wxy(1)=-wwz(1)*axisSign;     wxz(1)=wwy(1)*axisSign;      dxy(1)=-c(1)*axisSign;      dxz(1)=b(1)*axisSign;
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);      wxy(2)=-wwz(2)*axisSign;     wxz(2)=wwy(2)*axisSign;      dxy(2)=-c(2)*axisSign;      dxz(2)=b(2)*axisSign;
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);      wxy(3)=-wwz(3)*axisSign;     wxz(3)=wwy(3)*axisSign;      dxy(3)=-c(3)*axisSign;      dxz(3)=b(3)*axisSign;
               
        wwy=wwy.*edgeLength;    wwz=wwz.*edgeLength;    rw=rw.*edgeLength;
        for ii=1:3
            for jj=1:3
                Stt(ii,jj) = Stt(ii,jj) + Weights2D(kt) * (rw(ii) * imu * rw(jj));
                Szz(ii,jj) = Szz(ii,jj) + Weights2D(kt) * (dxy(ii) * imu * dxy(jj) + dxz(ii) * imu * dxz(jj));
                Ltz(ii,jj) = Ltz(ii,jj) + Weights2D(kt) * (wxy(ii) * imu * dxy(jj) + wxz(ii) * imu * dxz(jj));
                Lzt(ii,jj) = Lzt(ii,jj) + Weights2D(kt) * (dxy(ii) * imu * wxy(jj) + dxz(ii) * imu * wxz(jj));
                Ttt(ii,jj) = Ttt(ii,jj) + Weights2D(kt) * (wwy(ii) * epsilon *wwy(jj) + wwz(ii) * epsilon * wwz(jj));
                Tzz(ii,jj) = Tzz(ii,jj) + Weights2D(kt) * (simp(ii) * epsilon *simp(jj));
                Ktt(ii,jj) = Ktt(ii,jj) + Weights2D(kt) * (wxy(ii) * imu * wxy(jj) + wxz(ii) * imu *wxz(jj));
            end
        end
    end
    Stt=Stt*Ae;Szz=Szz*Ae;Ltz=Ltz*Ae;Lzt=Lzt*Ae;Ttt=Ttt*Ae;Tzz=Tzz*Ae;Ktt=Ktt*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=si*sj*Stt(ii,jj)-(omega^2)*si*sj*Ttt(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-si*sj*Ktt(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*Ltz(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counterM=counterM+1;IM(counterM)=vi_E;JM(counterM)=vj_E;M(counterM)=-Szz(ii,jj)+(omega^2)*Tzz(ii,jj);end
            if(vi_E~=0 && ej_E~=0),counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*Lzt(ii,jj);end
        end
    end
end
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM] = IVWE_Assembly_yBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==15)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 

                medium = element.Medium2D; epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;
    elseif(nargin==16)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 
                                                                                                                        freqIndex= varargin{16};
                medium = element.Medium2D; 
                if(medium.IsDispersive),epsilon = medium.Epsilon(freqIndex); mu = medium.Mu(freqIndex); waveImpedance = medium.WaveImpedance(freqIndex);
                else,epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;
                end
    end
    imu=(m0^-1)*mu^-1;epsilon=e0*epsilon;omega=2*pi*frequency;
    vertices =[toolboxModel.Vertices(element.Vertices)];
    edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];z=[vertices.Z];De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;    c(1)=(x(3)-x(2))/De;
    b(2)=(z(3)-z(1))/De;    c(2)=(x(1)-x(3))/De;
    b(3)=(z(1)-z(2))/De;    c(3)=(x(2)-x(1))/De;  
    %----------------------------------------------------------------------
    Stt=zeros(3,3);         Ttt=zeros(3,3);         Ltz=zeros(3,3);         Ktt=zeros(3,3);
    Szz=zeros(3,3);         Tzz=zeros(3,3);         Lzt=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);       wwz(1)=simp(1)*c(2)-simp(2)*c(1);
        wwx(2)=simp(2)*b(3)-simp(3)*b(2);       wwz(2)=simp(2)*c(3)-simp(3)*c(2);
        wwx(3)=simp(3)*b(1)-simp(1)*b(3);       wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);      wyx(1)=wwz(1)*axisSign;     wyz(1)=-wwx(1)*axisSign;      dyx(1)=c(1)*axisSign;      dyz(1)=-b(1)*axisSign;
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);      wyx(2)=wwz(2)*axisSign;     wyz(2)=-wwx(2)*axisSign;      dyx(2)=c(2)*axisSign;      dyz(2)=-b(2)*axisSign;
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);      wyx(3)=wwz(3)*axisSign;     wyz(3)=-wwx(3)*axisSign;      dyx(3)=c(3)*axisSign;      dyz(3)=-b(3)*axisSign;
               
        wwx=wwx.*edgeLength;    wwz=wwz.*edgeLength;    rw=rw.*edgeLength;
        for ii=1:3
            for jj=1:3
                Stt(ii,jj) = Stt(ii,jj) + Weights2D(kt) * (rw(ii) * imu * rw(jj));
                Szz(ii,jj) = Szz(ii,jj) + Weights2D(kt) * (dyx(ii) * imu * dyx(jj) + dyz(ii) * imu * dyz(jj));
                Ltz(ii,jj) = Ltz(ii,jj) + Weights2D(kt) * (wyx(ii) * imu * dyx(jj) + wyz(ii) * imu * dyz(jj));
                Lzt(ii,jj) = Lzt(ii,jj) + Weights2D(kt) * (dyx(ii) * imu * wyx(jj) + dyz(ii) * imu * wyz(jj));
                Ttt(ii,jj) = Ttt(ii,jj) + Weights2D(kt) * (wwx(ii) * epsilon *wwx(jj) + wwz(ii) * epsilon * wwz(jj));
                Tzz(ii,jj) = Tzz(ii,jj) + Weights2D(kt) * (simp(ii) * epsilon *simp(jj));
                Ktt(ii,jj) = Ktt(ii,jj) + Weights2D(kt) * (wyx(ii) * imu * wyx(jj) + wyz(ii) * imu *wyz(jj));
            end
        end
    end
    Stt=Stt*Ae;Szz=Szz*Ae;Ltz=Ltz*Ae;Lzt=Lzt*Ae;Ttt=Ttt*Ae;Tzz=Tzz*Ae;Ktt=Ktt*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=si*sj*Stt(ii,jj)-(omega^2)*si*sj*Ttt(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-si*sj*Ktt(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*Ltz(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counterM=counterM+1;IM(counterM)=vi_E;JM(counterM)=vj_E;M(counterM)=-Szz(ii,jj)+(omega^2)*Tzz(ii,jj);end
            if(vi_E~=0 && ej_E~=0),counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*Lzt(ii,jj);end
        end
    end
end
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM] = IVWE_Assembly_zBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==15)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 

                medium = element.Medium2D; epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;
    elseif(nargin==16)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 
                                                                                                                        freqIndex= varargin{16};
                medium = element.Medium2D; 
                if(medium.IsDispersive),epsilon = medium.Epsilon(freqIndex); mu = medium.Mu(freqIndex); waveImpedance = medium.WaveImpedance(freqIndex);
                else,epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;
                end
    end
    imu=(m0^-1)*mu^-1;epsilon=e0*epsilon;omega=2*pi*frequency;
    vertices =[toolboxModel.Vertices(element.Vertices)];
    edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];y=[vertices.Y];De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');Ae=abs(De)/2;
    b(1)=(y(2)-y(3))/De;    c(1)=(x(3)-x(2))/De;
    b(2)=(y(3)-y(1))/De;    c(2)=(x(1)-x(3))/De;
    b(3)=(y(1)-y(2))/De;    c(3)=(x(2)-x(1))/De;  
    %----------------------------------------------------------------------
    Stt=zeros(3,3);         Ttt=zeros(3,3);         Ltz=zeros(3,3);         Ktt=zeros(3,3);
    Szz=zeros(3,3);         Tzz=zeros(3,3);         Lzt=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);       wwy(1)=simp(1)*c(2)-simp(2)*c(1);
        wwx(2)=simp(2)*b(3)-simp(3)*b(2);       wwy(2)=simp(2)*c(3)-simp(3)*c(2);
        wwx(3)=simp(3)*b(1)-simp(1)*b(3);       wwy(3)=simp(3)*c(1)-simp(1)*c(3);
        
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);      wzx(1)=-wwy(1)*axisSign;     wzy(1)=wwx(1)*axisSign;      dzx(1)=-c(1)*axisSign;      dxy(1)=b(1)*axisSign;
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);      wzx(2)=-wwy(2)*axisSign;     wzy(2)=wwx(2)*axisSign;      dzx(2)=-c(2)*axisSign;      dxy(2)=b(2)*axisSign;
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);      wzx(3)=-wwy(3)*axisSign;     wzy(3)=wwx(3)*axisSign;      dzx(3)=-c(3)*axisSign;      dxy(3)=b(3)*axisSign;
               
        wwx=wwx.*edgeLength;    wwy=wwy.*edgeLength;    rw=rw.*edgeLength;
        for ii=1:3
            for jj=1:3
                Stt(ii,jj) = Stt(ii,jj) + Weights2D(kt) * (rw(ii) * imu * rw(jj));
                Szz(ii,jj) = Szz(ii,jj) + Weights2D(kt) * (dzx(ii) * imu * dzx(jj) + dzy(ii) * imu * dzy(jj));
                Ltz(ii,jj) = Ltz(ii,jj) + Weights2D(kt) * (wzx(ii) * imu * dzx(jj) + wzy(ii) * imu * dzy(jj));
                Lzt(ii,jj) = Lzt(ii,jj) + Weights2D(kt) * (dzx(ii) * imu * wzx(jj) + dzy(ii) * imu * wzy(jj));
                Ttt(ii,jj) = Ttt(ii,jj) + Weights2D(kt) * (wwx(ii) * epsilon *wwx(jj) + wwy(ii) * epsilon * wwy(jj));
                Tzz(ii,jj) = Tzz(ii,jj) + Weights2D(kt) * (simp(ii) * epsilon *simp(jj));
                Ktt(ii,jj) = Ktt(ii,jj) + Weights2D(kt) * (wzx(ii) * imu * wzx(jj) + wzy(ii) * imu *wzy(jj));
            end
        end
    end
    Stt=Stt*Ae;Szz=Szz*Ae;Ltz=Ltz*Ae;Lzt=Lzt*Ae;Ttt=Ttt*Ae;Tzz=Tzz*Ae;Ktt=Ktt*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=si*sj*Stt(ii,jj)-(omega^2)*si*sj*Ttt(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-si*sj*Ktt(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*Ltz(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counterM=counterM+1;IM(counterM)=vi_E;JM(counterM)=vj_E;M(counterM)=-Szz(ii,jj)+(omega^2)*Tzz(ii,jj);end
            if(vi_E~=0 && ej_E~=0),counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*Lzt(ii,jj);end
        end
    end
end
%==================== Anisotropic Assemblies ==============================
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM] = AVWE_Assembly_xBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==15)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 

                medium = element.Medium2D; epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;
    elseif(nargin==16)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 
                                                                                                                        freqIndex= varargin{16};
                medium = element.Medium2D; 
                if(medium.IsDispersive),epsilon = medium.Epsilon{freqIndex}; mu = medium.Mu{freqIndex}; waveImpedance = medium.WaveImpedance(freqIndex);
                else,epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;
                end
    end
    imu=(m0^-1)*mu^-1;epsilon=e0*epsilon;omega=2*pi*frequency;
    vertices =[toolboxModel.Vertices(element.Vertices)];
    edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];z=[vertices.Z];De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;    c(1)=(y(3)-y(2))/De;
    b(2)=(z(3)-z(1))/De;    c(2)=(y(1)-y(3))/De;
    b(3)=(z(1)-z(2))/De;    c(3)=(y(2)-y(1))/De;  
    %----------------------------------------------------------------------
    Stt=zeros(3,3);         Ttt=zeros(3,3);         Ltz=zeros(3,3);         Ktt=zeros(3,3);
    Szz=zeros(3,3);         Tzz=zeros(3,3);         Lzt=zeros(3,3);
    Stz=zeros(3,3);         Ttz=zeros(3,3);         LtA=zeros(3,3);
    Szt=zeros(3,3);         Tzt=zeros(3,3);         LtB=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwy(1)=simp(1)*b(2)-simp(2)*b(1);       wwz(1)=simp(1)*c(2)-simp(2)*c(1);
        wwy(2)=simp(2)*b(3)-simp(3)*b(2);       wwz(2)=simp(2)*c(3)-simp(3)*c(2);
        wwy(3)=simp(3)*b(1)-simp(1)*b(3);       wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);      wxy(1)=-wwz(1)*axisSign;     wxz(1)=wwy(1)*axisSign;      dxy(1)=-c(1)*axisSign;      dxz(1)=b(1)*axisSign;
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);      wxy(2)=-wwz(2)*axisSign;     wxz(2)=wwy(2)*axisSign;      dxy(2)=-c(2)*axisSign;      dxz(2)=b(2)*axisSign;
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);      wxy(3)=-wwz(3)*axisSign;     wxz(3)=wwy(3)*axisSign;      dxy(3)=-c(3)*axisSign;      dxz(3)=b(3)*axisSign;
               
        wwy=wwy.*edgeLength;    wwz=wwz.*edgeLength;    rw=rw.*edgeLength;
        for ii=1:3
            for jj=1:3
                Stt(ii,jj) = Stt(ii,jj) + Weights2D(kt) * (rw(ii) * imu(1,1) * rw(jj));
                Stz(ii,jj) = Stz(ii,jj) + Weights2D(kt) * (rw(ii) * imu(1,2) * dxy(jj) + rw(jj) * imu(1,3) * dxz(jj));
                Szt(ii,jj) = Szt(ii,jj) + Weights2D(kt) * (dxy(ii) * imu(2,1) * rw(jj) + dxz(jj) * imu(3,1) * rw(jj));
                Szz(ii,jj) = Szz(ii,jj) + Weights2D(kt) * (dxy(ii) * imu(2,2) * dxy(jj) + dxy(ii) * imu(2,3) *dxz(jj) + dxz(ii) * imu(3,2) * dxy(jj)  + dxz(ii) * imu(3,3) * dxz(jj));
                Ltz(ii,jj) = Ltz(ii,jj) + Weights2D(kt) * (wxy(ii) * imu(2,2) * dxy(jj) + wxy(ii) * imu(2,3) * dxz(jj) + wxz(ii) * imu(3,2) * dxy(jj) + wxz(ii) * imu(3,3) * dxz(jj));
                Lzt(ii,jj) = Lzt(ii,jj) + Weights2D(kt) * (dxy(ii) * imu(2,2) * wxy(jj) + dxy(ii) * imu(2,3) * wxz(jj) + dxz(ii) * imu(3,2) * wxy(jj) + dxz(ii) * imu(3,3) * wxz(jj));
                LtA(ii,jj) = LtA(ii,jj) + Weights2D(kt) * (rw(ii) * imu(1,2) * wxy(jj) + rw(ii) * imu(1,3) *wxz(jj));
                LtB(ii,jj) = LtB(ii,jj) + Weights2D(kt) * (wxy(ii) * imu(2,1) * rw(jj) + wxz(ii) * imu(3,1) * rw(jj));
                Ttt(ii,jj) = Ttt(ii,jj) + Weights2D(kt) * (wwy(ii) * epsilon(2,2) *wwy(jj) +wwy(ii) * epsilon(2,3) * wwz(jj) + wwz(ii) * epsilon(3,2) * wwy(jj) + wwz(ii) * epsilon(3,3) * wwz(jj));
                Tzz(ii,jj) = Tzz(ii,jj) + Weights2D(kt) * (simp(ii) * epsilon(1,1) *simp(jj));
                Ttz(ii,jj) = Ttz(ii,jj) + Weights2D(kt) * (wwy(ii) * epsilon(2,1) *simp(jj) +wwz(ii) * epsilon(3,1) * simp(jj));
                Tzt(ii,jj) = Tzt(ii,jj) + Weights2D(kt) * (simp(ii) * epsilon(1,2) *wwy(jj) +simp(ii) * epsilon(1,3) * wwz(jj));
                Ktt(ii,jj) = Ktt(ii,jj) + Weights2D(kt) * (wxy(ii) * imu(2,2) * wxy(jj) + wxy(ii) * imu(2,3) * wxz(jj) + wxz(ii) * imu(3,3) *wxz(jj) +wxz(ii) * imu(3,2) * wxy(jj));
            end
        end
    end
    Stt=Stt*Ae;     Ttt=Ttt*Ae;     Ltz=Ltz*Ae;         Ktt=Ktt*Ae;
    Szz=Szz*Ae;     Tzz=Tzz*Ae;     Lzt=Lzt*Ae;
    Stz=Stz*Ae;     Ttz=Ttz*Ae;     LtA=LtA*Ae;
    Szt=Szt*Ae;     Tzt=Tzt*Ae;     LtB=LtB*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=si*sj*Stt(ii,jj)-(omega^2)*si*sj*Ttt(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-si*sj*Ktt(ii,jj);
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=ej_E;L(counterL)=si*sj*LtB(ii,jj)-si*sj*LtA(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*Ltz(ii,jj);
                                   counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=vj_E;M(counterM)=si*Stz(ii,jj) - (omega^2)*si*Ttz(ii,jj);
            end
            if(vi_E~=0 && vj_E~=0),counterM=counterM+1;IM(counterM)=vi_E;JM(counterM)=vj_E;M(counterM)=-Szz(ii,jj)+(omega^2)*Tzz(ii,jj);end
            if(vi_E~=0 && ej_E~=0),counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*Lzt(ii,jj);
                                   counterM=counterM+1;IM(counterM)=vi_E;JM(counterM)=ej_E;M(counterM)=-sj*Szt(ii,jj) + (omega^2)*sj*Tzt(ii,jj);
            end
        end
    end
end
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM] = AVWE_Assembly_yBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==15)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 

                medium = element.Medium2D; epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;
    elseif(nargin==16)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 
                                                                                                                        freqIndex= varargin{16};
                medium = element.Medium2D; 
                if(medium.IsDispersive),epsilon = medium.Epsilon{freqIndex}; mu = medium.Mu{freqIndex}; waveImpedance = medium.WaveImpedance(freqIndex);
                else,epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;
                end
    end
    imu=(m0^-1)*mu^-1;epsilon=e0*epsilon;omega=2*pi*frequency;
    vertices =[toolboxModel.Vertices(element.Vertices)];
    edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];z=[vertices.Z];De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;    c(1)=(x(3)-x(2))/De;
    b(2)=(z(3)-z(1))/De;    c(2)=(x(1)-x(3))/De;
    b(3)=(z(1)-z(2))/De;    c(3)=(x(2)-x(1))/De;  
    %----------------------------------------------------------------------
    Stt=zeros(3,3);         Ttt=zeros(3,3);         Ltz=zeros(3,3);         Ktt=zeros(3,3);
    Szz=zeros(3,3);         Tzz=zeros(3,3);         Lzt=zeros(3,3);
    Stz=zeros(3,3);         Ttz=zeros(3,3);         LtA=zeros(3,3);
    Szt=zeros(3,3);         Tzt=zeros(3,3);         LtB=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);       wwz(1)=simp(1)*c(2)-simp(2)*c(1);
        wwx(2)=simp(2)*b(3)-simp(3)*b(2);       wwz(2)=simp(2)*c(3)-simp(3)*c(2);
        wwx(3)=simp(3)*b(1)-simp(1)*b(3);       wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);      wyx(1)=wwz(1)*axisSign;     wyz(1)=-wwx(1)*axisSign;      dyx(1)=c(1)*axisSign;      dyz(1)=-b(1)*axisSign;
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);      wyx(2)=wwz(2)*axisSign;     wyz(2)=-wwx(2)*axisSign;      dyx(2)=c(2)*axisSign;      dyz(2)=-b(2)*axisSign;
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);      wyx(3)=wwz(3)*axisSign;     wyz(3)=-wwx(3)*axisSign;      dyx(3)=c(3)*axisSign;      dyz(3)=-b(3)*axisSign;
               
        wwx=wwx.*edgeLength;    wwz=wwz.*edgeLength;    rw=rw.*edgeLength;
        for ii=1:3
            for jj=1:3
                Stt(ii,jj) = Stt(ii,jj) + Weights2D(kt) * (rw(ii) * imu(2,2) * rw(jj));
                Stz(ii,jj) = Stz(ii,jj) + Weights2D(kt) * (rw(ii) * imu(2,1) * dyx(jj) + rw(jj) * imu(2,3) * dyz(jj));
                Szt(ii,jj) = Szt(ii,jj) + Weights2D(kt) * (dyx(ii) * imu(1,2) * rw(jj) + dyz(jj) * imu(3,2) * rw(jj));
                Szz(ii,jj) = Szz(ii,jj) + Weights2D(kt) * (dyx(ii) * imu(1,1) * dyx(jj) + dyx(ii) * imu(1,3) *dyz(jj) + dyz(ii) * imu(3,1) * dyx(jj)  + dyz(ii) * imu(3,3) * dyz(jj));
                Ltz(ii,jj) = Ltz(ii,jj) + Weights2D(kt) * (wyx(ii) * imu(1,1) * dyx(jj) + wyx(ii) * imu(1,3) * dyz(jj) + wyz(ii) * imu(3,1) * dyx(jj) + wyz(ii) * imu(3,3) * dyz(jj));
                Lzt(ii,jj) = Lzt(ii,jj) + Weights2D(kt) * (dyx(ii) * imu(1,1) * wyx(jj) + dyx(ii) * imu(1,3) * wyz(jj) + dyz(ii) * imu(3,1) * wyx(jj) + dyz(ii) * imu(3,3) * wyz(jj));
                LtA(ii,jj) = LtA(ii,jj) + Weights2D(kt) * (rw(ii) * imu(2,1) * wyx(jj) + rw(ii) * imu(2,3) *wyz(jj));
                LtB(ii,jj) = LtB(ii,jj) + Weights2D(kt) * (wyx(ii) * imu(1,2) * rw(jj) + wyz(ii) * imu(3,2) * rw(jj));
                Ttt(ii,jj) = Ttt(ii,jj) + Weights2D(kt) * (wwx(ii) * epsilon(1,1) *wwx(jj) +wwx(ii) * epsilon(1,3) * wwz(jj) + wwz(ii) * epsilon(3,1) * wwx(jj) + wwz(ii) * epsilon(3,3) * wwz(jj));
                Tzz(ii,jj) = Tzz(ii,jj) + Weights2D(kt) * (simp(ii) * epsilon(2,2) *simp(jj));
                Ttz(ii,jj) = Ttz(ii,jj) + Weights2D(kt) * (wwx(ii) * epsilon(1,2) *simp(jj) +wwz(ii) * epsilon(3,2) * simp(jj));
                Tzt(ii,jj) = Tzt(ii,jj) + Weights2D(kt) * (simp(ii) * epsilon(2,1) *wwx(jj) +simp(ii) * epsilon(2,3) * wwz(jj));
                Ktt(ii,jj) = Ktt(ii,jj) + Weights2D(kt) * (wyx(ii) * imu(1,1) * wyx(jj) + wyx(ii) * imu(1,3) * wyz(jj) + wyz(ii) * imu(3,3) *wyz(jj) +wyz(ii) * imu(3,1) * wyx(jj));
            end
        end
    end
    Stt=Stt*Ae;     Ttt=Ttt*Ae;     Ltz=Ltz*Ae;         Ktt=Ktt*Ae;
    Szz=Szz*Ae;     Tzz=Tzz*Ae;     Lzt=Lzt*Ae;
    Stz=Stz*Ae;     Ttz=Ttz*Ae;     LtA=LtA*Ae;
    Szt=Szt*Ae;     Tzt=Tzt*Ae;     LtB=LtB*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=si*sj*Stt(ii,jj)-(omega^2)*si*sj*Ttt(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-si*sj*Ktt(ii,jj);
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=ej_E;L(counterL)=si*sj*LtB(ii,jj)-si*sj*LtA(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*Ltz(ii,jj);
                                   counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=vj_E;M(counterM)=si*Stz(ii,jj) - (omega^2)*si*Ttz(ii,jj);
            end
            if(vi_E~=0 && vj_E~=0),counterM=counterM+1;IM(counterM)=vi_E;JM(counterM)=vj_E;M(counterM)=-Szz(ii,jj)+(omega^2)*Tzz(ii,jj);end
            if(vi_E~=0 && ej_E~=0),counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*Lzt(ii,jj);
                                   counterM=counterM+1;IM(counterM)=vi_E;JM(counterM)=ej_E;M(counterM)=-sj*Szt(ii,jj) + (omega^2)*sj*Tzt(ii,jj);
            end
        end
    end
end
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM] = AVWE_Assembly_zBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==15)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 

                medium = element.Medium2D; epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;
    elseif(nargin==16)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 
                                                                                                                        freqIndex= varargin{16};
                medium = element.Medium2D; 
                if(medium.IsDispersive),epsilon = medium.Epsilon{freqIndex}; mu = medium.Mu{freqIndex}; waveImpedance = medium.WaveImpedance(freqIndex);
                else,epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;
                end
    end
    imu=(m0^-1)*mu^-1;epsilon=e0*epsilon;omega=2*pi*frequency;
    vertices =[toolboxModel.Vertices(element.Vertices)];
    edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];y=[vertices.Y];De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');Ae=abs(De)/2;
    b(1)=(y(2)-y(3))/De;    c(1)=(x(3)-x(2))/De;
    b(2)=(y(3)-y(1))/De;    c(2)=(x(1)-x(3))/De;
    b(3)=(y(1)-y(2))/De;    c(3)=(x(2)-x(1))/De;  
    %----------------------------------------------------------------------
    Stt=zeros(3,3);         Ttt=zeros(3,3);         Ltz=zeros(3,3);         Ktt=zeros(3,3);
    Szz=zeros(3,3);         Tzz=zeros(3,3);         Lzt=zeros(3,3);
    Stz=zeros(3,3);         Ttz=zeros(3,3);         LtA=zeros(3,3);
    Szt=zeros(3,3);         Tzt=zeros(3,3);         LtB=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);       wwy(1)=simp(1)*c(2)-simp(2)*c(1);
        wwx(2)=simp(2)*b(3)-simp(3)*b(2);       wwy(2)=simp(2)*c(3)-simp(3)*c(2);
        wwx(3)=simp(3)*b(1)-simp(1)*b(3);       wwy(3)=simp(3)*c(1)-simp(1)*c(3);
        
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);      wzx(1)=-wwy(1)*axisSign;     wzy(1)=wwx(1)*axisSign;      dzx(1)=-c(1)*axisSign;      dzy(1)=b(1)*axisSign;
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);      wzx(2)=-wwy(2)*axisSign;     wzy(2)=wwx(2)*axisSign;      dzx(2)=-c(2)*axisSign;      dzy(2)=b(2)*axisSign;
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);      wzx(3)=-wwy(3)*axisSign;     wzy(3)=wwx(3)*axisSign;      dzx(3)=-c(3)*axisSign;      dzy(3)=b(3)*axisSign;
               
        wwx=wwx.*edgeLength;    wwy=wwy.*edgeLength;    rw=rw.*edgeLength;
        for ii=1:3
            for jj=1:3
                Stt(ii,jj) = Stt(ii,jj) + Weights2D(kt) * (rw(ii) * imu(3,3) * rw(jj));
                Stz(ii,jj) = Stz(ii,jj) + Weights2D(kt) * (rw(ii) * imu(3,1) * dzx(jj) + rw(jj) * imu(3,2) * dzy(jj));
                Szt(ii,jj) = Szt(ii,jj) + Weights2D(kt) * (dzx(ii) * imu(1,3) * rw(jj) + dzy(jj) * imu(2,3) * rw(jj));
                Szz(ii,jj) = Szz(ii,jj) + Weights2D(kt) * (dzx(ii) * imu(1,1) * dzx(jj) + dzx(ii) * imu(1,2) *dzy(jj) + dzy(ii) * imu(2,1) * dzx(jj)  + dzy(ii) * imu(2,2) * dzy(jj));
                Ltz(ii,jj) = Ltz(ii,jj) + Weights2D(kt) * (wzx(ii) * imu(1,1) * dzx(jj) + wzx(ii) * imu(1,2) * dzy(jj) + wzy(ii) * imu(2,1) * dzx(jj) + wzy(ii) * imu(2,2) * dzy(jj));
                Lzt(ii,jj) = Lzt(ii,jj) + Weights2D(kt) * (dzx(ii) * imu(1,1) * wzx(jj) + dzx(ii) * imu(1,2) * wzy(jj) + dzy(ii) * imu(2,1) * wzx(jj) + dzy(ii) * imu(2,2) * wzy(jj));
                LtA(ii,jj) = LtA(ii,jj) + Weights2D(kt) * (rw(ii) * imu(3,1) * wzx(jj) + rw(ii) * imu(3,2) *wzy(jj));
                LtB(ii,jj) = LtB(ii,jj) + Weights2D(kt) * (wzx(ii) * imu(1,3) * rw(jj) + wzy(ii) * imu(2,3) * rw(jj));
                Ttt(ii,jj) = Ttt(ii,jj) + Weights2D(kt) * (wwx(ii) * epsilon(1,1) *wwx(jj) +wwx(ii) * epsilon(1,2) * wwy(jj) + wwy(ii) * epsilon(2,1) * wwx(jj) + wwy(ii) * epsilon(2,2) * wwy(jj));
                Tzz(ii,jj) = Tzz(ii,jj) + Weights2D(kt) * (simp(ii) * epsilon(3,3) *simp(jj));
                Ttz(ii,jj) = Ttz(ii,jj) + Weights2D(kt) * (wwx(ii) * epsilon(1,3) *simp(jj) +wwy(ii) * epsilon(2,3) * simp(jj));
                Tzt(ii,jj) = Tzt(ii,jj) + Weights2D(kt) * (simp(ii) * epsilon(3,1) *wwx(jj) +simp(ii) * epsilon(3,2) * wwy(jj));
                Ktt(ii,jj) = Ktt(ii,jj) + Weights2D(kt) * (wzx(ii) * imu(1,1) * wzx(jj) + wzx(ii) * imu(1,2) * wzy(jj) + wzy(ii) * imu(2,2) *wzy(jj) +wzy(ii) * imu(2,1) * wxy(jj));
            end
        end
    end
    Stt=Stt*Ae;     Ttt=Ttt*Ae;     Ltz=Ltz*Ae;         Ktt=Ktt*Ae;
    Szz=Szz*Ae;     Tzz=Tzz*Ae;     Lzt=Lzt*Ae;
    Stz=Stz*Ae;     Ttz=Ttz*Ae;     LtA=LtA*Ae;
    Szt=Szt*Ae;     Tzt=Tzt*Ae;     LtB=LtB*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=si*sj*Stt(ii,jj)-(omega^2)*si*sj*Ttt(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-si*sj*Ktt(ii,jj);
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=ej_E;L(counterL)=si*sj*LtB(ii,jj)-si*sj*LtA(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*Ltz(ii,jj);
                                   counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=vj_E;M(counterM)=si*Stz(ii,jj) - (omega^2)*si*Ttz(ii,jj);
            end
            if(vi_E~=0 && vj_E~=0),counterM=counterM+1;IM(counterM)=vi_E;JM(counterM)=vj_E;M(counterM)=-Szz(ii,jj)+(omega^2)*Tzz(ii,jj);end
            if(vi_E~=0 && ej_E~=0),counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*Lzt(ii,jj);
                                   counterM=counterM+1;IM(counterM)=vi_E;JM(counterM)=ej_E;M(counterM)=-sj*Szt(ii,jj) + (omega^2)*sj*Tzt(ii,jj);
            end
        end
    end
end
%=================== Bianisotropic Assemblied =============================
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM] = BVWE_Assembly_xBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==15)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 

                medium = element.Medium2D; epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;
                ksi = medium.Ksi; zeta = medium.Zita; 
    elseif(nargin==16)
                IK = varargin{1};       IL = varargin{4};       IM = varargin{7};       counterK = varargin{10};        toolboxModel = varargin{13};       
                JK = varargin{2};       JL = varargin{5};       JM = varargin{8};       counterL = varargin{11};        element = varargin{14};
                K  = varargin{3};       L  = varargin{6};       M  = varargin{9};       counterM = varargin{12};        axisSign = varargin{15}; 
                                                                                                                        freqIndex= varargin{16};
                medium = element.Medium2D; 
                if(medium.IsDispersive),epsilon = medium.Epsilon{freqIndex}; mu = medium.Mu{freqIndex}; waveImpedance = medium.WaveImpedance(freqIndex);ksi=medium.Ksi{freqIndex};zeta = medium.Zita{freqIndex};
                else,epsilon = medium.Epsilon; mu = medium.Mu; waveImpedance = medium.WaveImpedance;ksi=medium.Ksi;zeta=medium.Zita;
                end
    end
    imu=(m0^-1)*mu^-1;epsilon=e0*epsilon;omega=2*pi*frequency;phi=imu*zeta;theta=epsilon-ksi*imu*zeta;psi=ksi*imu;
    vertices =[toolboxModel.Vertices(element.Vertices)];
    edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];z=[vertices.Z];De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;    c(1)=(y(3)-y(2))/De;
    b(2)=(z(3)-z(1))/De;    c(2)=(y(1)-y(3))/De;
    b(3)=(z(1)-z(2))/De;    c(3)=(y(2)-y(1))/De;  
    %----------------------------------------------------------------------
    Stt=zeros(3,3);         Ttt=zeros(3,3);         LtzA=zeros(3,3);     LtA=zeros(3,3);     MttA=zeros(3,3);       MtzA=zeros(3,3);      Ktt=zeros(3,3);
    Szz=zeros(3,3);         Tzz=zeros(3,3);         LztA=zeros(3,3);     LtB=zeros(3,3);     MttB=zeros(3,3);       MtzB=zeros(3,3);
    Stz=zeros(3,3);         Ttz=zeros(3,3);         LtzB=zeros(3,3);     LtC=zeros(3,3);     MzzA=zeros(3,3);       MztA=zeros(3,3);   
    Szt=zeros(3,3);         Tzt=zeros(3,3);         LztB=zeros(3,3);     LtD=zeros(3,3);     MzzB=zeros(3,3);       MztB=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwy(1)=simp(1)*b(2)-simp(2)*b(1);       wwz(1)=simp(1)*c(2)-simp(2)*c(1);
        wwy(2)=simp(2)*b(3)-simp(3)*b(2);       wwz(2)=simp(2)*c(3)-simp(3)*c(2);
        wwy(3)=simp(3)*b(1)-simp(1)*b(3);       wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);      wxy(1)=-wwz(1)*axisSign;     wxz(1)=wwy(1)*axisSign;      dxy(1)=-c(1)*axisSign;      dxz(1)=b(1)*axisSign;
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);      wxy(2)=-wwz(2)*axisSign;     wxz(2)=wwy(2)*axisSign;      dxy(2)=-c(2)*axisSign;      dxz(2)=b(2)*axisSign;
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);      wxy(3)=-wwz(3)*axisSign;     wxz(3)=wwy(3)*axisSign;      dxy(3)=-c(3)*axisSign;      dxz(3)=b(3)*axisSign;
               
        wwy=wwy.*edgeLength;    wwz=wwz.*edgeLength;    rw=rw.*edgeLength;
        for ii=1:3
            for jj=1:3
                Stt(ii,jj) = Stt(ii,jj) + Weights2D(kt) * (rw(ii) * imu(1,1) * rw(jj));
                Stz(ii,jj) = Stz(ii,jj) + Weights2D(kt) * (rw(ii) * imu(1,2) * dxy(jj) + rw(jj) * imu(1,3) * dxz(jj));
                Szt(ii,jj) = Szt(ii,jj) + Weights2D(kt) * (dxy(ii) * imu(2,1) * rw(jj) + dxz(jj) * imu(3,1) * rw(jj));
                Szz(ii,jj) = Szz(ii,jj) + Weights2D(kt) * (dxy(ii) * imu(2,2) * dxy(jj) + dxy(ii) * imu(2,3) *dxz(jj) + dxz(ii) * imu(3,2) * dxy(jj)  + dxz(ii) * imu(3,3) * dxz(jj));
                LtzA(ii,jj) = LtzA(ii,jj) + Weights2D(kt) * (wxy(ii) * imu(2,2) * dxy(jj) + wxy(ii) * imu(2,3) * dxz(jj) + wxz(ii) * imu(3,2) * dxy(jj) + wxz(ii) * imu(3,3) * dxz(jj));
                LztA(ii,jj) = LztA(ii,jj) + Weights2D(kt) * (dxy(ii) * imu(2,2) * wxy(jj) + dxy(ii) * imu(2,3) * wxz(jj) + dxz(ii) * imu(3,2) * wxy(jj) + dxz(ii) * imu(3,3) * wxz(jj));
                LtA(ii,jj) = LtA(ii,jj) + Weights2D(kt) * (rw(ii) * imu(1,2) * wxy(jj) + rw(ii) * imu(1,3) *wxz(jj));
                LtB(ii,jj) = LtB(ii,jj) + Weights2D(kt) * (wxy(ii) * imu(2,1) * rw(jj) + wxz(ii) * imu(3,1) * rw(jj));
                Ttt(ii,jj) = Ttt(ii,jj) + Weights2D(kt) * (wwy(ii) * epsilon(2,2) *wwy(jj) +wwy(ii) * epsilon(2,3) * wwz(jj) + wwz(ii) * epsilon(3,2) * wwy(jj) + wwz(ii) * epsilon(3,3) * wwz(jj));
                Tzz(ii,jj) = Tzz(ii,jj) + Weights2D(kt) * (simp(ii) * epsilon(1,1) *simp(jj));
                Ttz(ii,jj) = Ttz(ii,jj) + Weights2D(kt) * (wwy(ii) * epsilon(2,1) *simp(jj) +wwz(ii) * epsilon(3,1) * simp(jj));
                Tzt(ii,jj) = Tzt(ii,jj) + Weights2D(kt) * (simp(ii) * epsilon(1,2) *wwy(jj) +simp(ii) * epsilon(1,3) * wwz(jj));
                Ktt(ii,jj) = Ktt(ii,jj) + Weights2D(kt) * (wxy(ii) * imu(2,2) * wxy(jj) + wxy(ii) * imu(2,3) * wxz(jj) + wxz(ii) * imu(3,3) *wxz(jj) +wxz(ii) * imu(3,2) * wxy(jj));
                %----------------------------------------------------------
                LtC(ii,jj) = LtC(ii,jj) + Weights2D(kt) * (wxy(ii) * phi(2,2) *wwy(jj) + wxy(ii) * phi(2,3) * wwz(jj) + wxz(ii) * phi(3,2) * wwy(jj) + wxz(ii) * phi(3,3) * wwz(jj));
                LtD(ii,jj) = LtD(ii,jj) + Weights2D(kt) * ();
                LtzB(ii,jj) = LtzB(ii,jj) + Weights2D(kt)*();
                LztB(ii,jj) = LztB(ii,jj) + Weights2D(kt)*();
                MttA(ii,jj) = MttA(ii,jj) + Weights2D(kt)*();
                MttB(ii,jj) = MttB(ii,jj) + Weights2D(kt)*();
                MzzA(ii,jj) = MzzA(ii,jj) + Weights2D(kt)*();
                MzzB(ii,jj) = MzzB(ii,jj) + Weights2D(kt)*();
                MtzA(ii,jj) = MtzA(ii,jj) + Weights2D(kt)*();
                MtzB(ii,jj) = MtzB(ii,jj) + Weights2D(kt)*();
                MztA(ii,jj) = MztA(ii,jj) + Weights2D(kt)*();
                MztB(ii,jj) = MztB(ii,jj) + Weights2D(kt)*();
            end
        end
    end
    Stt=Stt*Ae;     Ttt=Ttt*Ae;     LtzA=LtzA*Ae;         Ktt=Ktt*Ae;
    Szz=Szz*Ae;     Tzz=Tzz*Ae;     LztA=LztA*Ae;
    Stz=Stz*Ae;     Ttz=Ttz*Ae;     LtA=LtA*Ae;
    Szt=Szt*Ae;     Tzt=Tzt*Ae;     LtB=LtB*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=si*sj*Stt(ii,jj)-(omega^2)*si*sj*Ttt(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-si*sj*Ktt(ii,jj);
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=ej_E;L(counterL)=si*sj*LtB(ii,jj)-si*sj*LtA(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*LtzA(ii,jj);
                                   counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=vj_E;M(counterM)=si*Stz(ii,jj) - (omega^2)*si*Ttz(ii,jj);
            end
            if(vi_E~=0 && vj_E~=0),counterM=counterM+1;IM(counterM)=vi_E;JM(counterM)=vj_E;M(counterM)=-Szz(ii,jj)+(omega^2)*Tzz(ii,jj);end
            if(vi_E~=0 && ej_E~=0),counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*LztA(ii,jj);
                                   counterM=counterM+1;IM(counterM)=vi_E;JM(counterM)=ej_E;M(counterM)=-sj*Szt(ii,jj) + (omega^2)*sj*Tzt(ii,jj);
            end
        end
    end
end