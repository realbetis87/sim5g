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
function [AssembledSystem] = VWE_Assembly(toolboxModel,AssembledSystem,boundaryIndices),frequency=toolboxModel.Frequency;
    if(frequency.NF==1)
        if(isvector(boundaryIndices)),boundary=toolboxModel.Boundaries(boundaryIndices(1));
            switch boundary.Axis
                case  1,[AssembledSystem] = VWE_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices);
                case -1,[AssembledSystem] = VWE_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices);
                case  2,[AssembledSystem] = VWE_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices);
                case -2,[AssembledSystem] = VWE_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices);
                case  3,[AssembledSystem] = VWE_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices);
                case -3,[AssembledSystem] = VWE_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices);
            end
        else,boundary=toolboxModel.Boundaries(boundaryIndices);
             switch boundary.Axis
                case  1,[AssembledSystem] = VWE_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices);
                case -1,[AssembledSystem] = VWE_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices);
                case  2,[AssembledSystem] = VWE_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices);
                case -2,[AssembledSystem] = VWE_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices);
                case  3,[AssembledSystem] = VWE_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices);
                case -3,[AssembledSystem] = VWE_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices);
            end
        end
    else
        if(isvector(boundaryIndices)),boundary=toolboxModel.Boundaries(boundaryIndices(1));
            switch boundary.Axis
                case 1, for ii=1:frequency.NF,[AssembledSystem] = VWE_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices,ii);end
                case -1,for ii=1:frequency.NF,[AssembledSystem] = VWE_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices,ii);end
                case 2, for ii=1:frequency.NF,[AssembledSystem] = VWE_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices,ii);end
                case -2,for ii=1:frequency.NF,[AssembledSystem] = VWE_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices,ii);end
                case 3, for ii=1:frequency.NF,[AssembledSystem] = VWE_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices,ii);end
                case -3,for ii=1:frequency.NF,[AssembledSystem] = VWE_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices,ii);end
            end
        else,boundary=toolboxModel.Boundaries(boundaryIndices);
             switch boundary.Axis
                 case 1, for ii=1:frequency.NF,[AssembledSystem] = VWE_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices,ii);end
                 case -1,for ii=1:frequency.NF,[AssembledSystem] = VWE_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices,ii);end
                 case 2, for ii=1:frequency.NF,[AssembledSystem] = VWE_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices,ii);end
                 case -2,for ii=1:frequency.NF,[AssembledSystem] = VWE_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices,ii);end
                 case 3, for ii=1:frequency.NF,[AssembledSystem] = VWE_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices,ii);end
                 case -3,for ii=1:frequency.NF,[AssembledSystem] = VWE_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices,ii);end
            end
        end
    end
end
%---------------------------- x Boundaries --------------------------------
function [AssembledSystem] = VWE_Assembly_xBoundary(varargin)
    if(nargin==3),toolboxModel=varargin{1};AssembledSystem=varargin{2};boundaryIndices=varargin{3};
            N=AssembledSystem.DimEt+AssembledSystem.DimEn;
            IK=zeros(120*9*N,1);JK=zeros(120*9*N,1);K=zeros(120*9*N,1);IL=zeros(120*9*N,1);JL=zeros(120*9*N,1);L=zeros(120*9*N,1);
            IM=zeros(120*9*N,1);JM=zeros(120*9*N,1);M=zeros(120*9*N,1);counterK=0;counterL=0;counterM=0;
            for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
                for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                    switch medium.Type
                        case "Iso", [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_xBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element);
                        case "Anis",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_xBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element);
                        case "Bian",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_xBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element);
                    end
                end
            end
            Nv=AssembledSystem.DimEn;
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

            AssembledSystem.Matrix_A=SA;    AssembledSystem.Matrix_B=SB;
    elseif(nargin==4),toolboxModel=varargin{1};AssembledSystem=varargin{2};boundaryIndices=varargin{3};FrequencyIndex=varargin{4};
        N=AssembledSystem.DimEt+AssembledSystem.DimEn;
            IK=zeros(120*9*N,1);JK=zeros(120*9*N,1);K=zeros(120*9*N,1);IL=zeros(120*9*N,1);JL=zeros(120*9*N,1);L=zeros(120*9*N,1);
            IM=zeros(120*9*N,1);JM=zeros(120*9*N,1);M=zeros(120*9*N,1);counterK=0;counterL=0;counterM=0;
            for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
                for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                    switch medium.Type
                        case "Iso", [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_xBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,FrequencyIndex);
                        case "Anis",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_xBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,FrequencyIndex);
                        case "Bian",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_xBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,element,FrequencyIndex);
                    end
                end
            end
            Nv=AssembledSystem.DimEn;
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

            AssembledSystem.Matrix_A{FrequencyIndex}=SA;    AssembledSystem.Matrix_B{FrequencyIndex}=SB;
    end
end
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_xBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==14),IK=varargin{1};      IL=varargin{4};     IM=varargin{7};         counterK=varargin{10};          toolboxModel=varargin{13};
                   JK=varargin{2};      JL=varargin{5};     JM=varargin{8};         counterL=varargin{11};          element=varargin{14};
                   K=varargin{3};       L=varargin{6};      M=varargin{9};          counterM=varargin{12};
                   
                   medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;waveImpedance=medium.WaveImpedance;
                   frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==15),IK=varargin{1};   IL=varargin{4};      IM=varargin{7};       counterK=varargin{10};      toolboxModel=varargin{13};
                       JK=varargin{2};   JL=varargin{5};      JM=varargin{8};       counterL=varargin{11};      element=varargin{14};
                       K=varargin{3};    L=varargin{6};       M=varargin{9};        counterM=varargin{12};      FrequencyIndex=varargin{15};
                  medium=element.Medium2D;
                   if(medium.IsDispersive),epsilon=medium.Epsilon{FrequencyIndex};mu=medium.Mu{FrequencyIndex};waveImpedance=medium.WaveImpedance{FrequencyIndex};
                   else,epsilon=medium.Epsilon;mu=medium.Mu;waveImpedance=medium.WaveImpedance;
                   end,frequency=toolboxModel.Frequency.Frequency(FrequencyIndex);
    end,im=(m0^-1)*mu^-1;epsilon=e0*epsilon;omega=2*pi*frequency;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];
    edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];z=[vertices.Z];De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;    c(1)=(y(3)-y(2))/De;
    b(2)=(z(3)-z(1))/De;    c(2)=(y(1)-y(3))/De;
    b(3)=(z(1)-z(2))/De;    c(3)=(y(2)-y(1))/De;   
    T11=zeros(3,3);T15=zeros(3,3);T16=zeros(3,3);T18=zeros(3,3);T19=zeros(3,3);
    T21=zeros(3,3);T23=zeros(3,3);Tsge=zeros(3,3);Tsgv=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwy(1)=simp(1)*b(2)-simp(2)*b(1);       wwz(1)=simp(1)*c(2)-simp(2)*c(1);
        wwy(2)=simp(2)*b(3)-simp(3)*b(2);       wwz(2)=simp(2)*c(3)-simp(3)*c(2);
        wwy(3)=simp(3)*b(1)-simp(1)*b(3);       wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);      wxy(1)=-wwz(1);     wxz(1)=wwy(1);      dxy(1)=-c(1);      dxz(1)=b(1);
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);      wxy(2)=-wwz(2);     wxz(2)=wwy(2);      dxy(2)=-c(2);      dxz(2)=b(2);
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);      wxy(3)=-wwz(3);     wxz(3)=wwy(3);      dxy(3)=-c(3);      dxz(3)=b(3);
               
        wwy=wwy.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        for ii=1:3
            for jj=1:3,T11(ii,jj)=T11(ii,jj)+Weights2D(kt)*(rw(ii)*im*rw(jj));
                       T15(ii,jj)=T15(ii,jj)+Weights2D(kt)*(wxy(ii)*im*wxy(jj)+wxz(ii)*im*wxz(jj));
                       T16(ii,jj)=T16(ii,jj)+Weights2D(kt)*(wxy(ii)*im*dxy(jj)+wxz(ii)*im*dxz(jj));
                       T18(ii,jj)=T18(ii,jj)+Weights2D(kt)*(dxy(ii)*im*wxy(jj)+dxz(ii)*im*wxz(jj));%
                       T19(ii,jj)=T19(ii,jj)+Weights2D(kt)*(dxy(ii)*im*dxy(jj)+dxz(ii)*im*dxz(jj));
                       T21(ii,jj)=T21(ii,jj)+Weights2D(kt)*(wwy(ii)*epsilon*wwy(jj)+wwz(ii)*epsilon*wwz(jj));
                       T23(ii,jj)=T23(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon*simp(jj));
            end
        end
    end
    T11=T11*Ae;T15=T15*Ae;T16=T16*Ae;T18=T18*Ae;T19=T19*Ae;T21=T21*Ae;T23=T23*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-(omega^2)*si*sj*T21(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=si*sj*T11(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=1i*k0*waveImpedance*si*sj*Tsge(ii,jj);
                                   %------------------------------------------------------------------------------------
                                   counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=-si*sj*T15(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*T16(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=(omega^2)*T23(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=-T19(ii,jj);
            end
            if(vi_E~=0 && ej_E~=0),counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*T18(ii,jj);end
        end
    end
end
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_xBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==14),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};
                   medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;
                   frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==15),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};FrequencyIndex=varargin{15};medium=element.Medium2D;
                   if(medium.IsDispersive),epsilon=medium.Epsilon{FrequencyIndex};mu=medium.Mu{FrequencyIndex};
                   else,epsilon=medium.Epsilon;mu=medium.Mu;
                   end,frequency=toolboxModel.Frequency.Frequency(FrequencyIndex);
    end,im=(m0^-1)*mu^-1;epsilon=e0*epsilon;omega=2*pi*frequency;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];z=[vertices.Z];De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;b(2)=(z(3)-z(1))/De;b(3)=(z(1)-z(2))/De;c(1)=(y(3)-y(2))/De;c(2)=(y(1)-y(3))/De;c(3)=(y(2)-y(1))/De;   
    T11=zeros(3,3);T12=zeros(3,3);T13=zeros(3,3);T14=zeros(3,3);T15=zeros(3,3);T16=zeros(3,3);T17=zeros(3,3);T18=zeros(3,3);T19=zeros(3,3);
    T21=zeros(3,3);T22=zeros(3,3);T23=zeros(3,3);T24=zeros(3,3);
    Tsge=zeros(3,3);Tsgv=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwy(1)=simp(1)*b(2)-simp(2)*b(1);wwy(2)=simp(2)*b(3)-simp(3)*b(2);wwy(3)=simp(3)*b(1)-simp(1)*b(3);
        wwz(1)=simp(1)*c(2)-simp(2)*c(1);wwz(2)=simp(2)*c(3)-simp(3)*c(2);wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wxy(1)=-wwz(1);wxz(1)=wwy(1);wxy(2)=-wwz(2);wxz(2)=wwy(2);wxy(3)=-wwz(3);wxz(3)=wwy(3);
        wwy=wwy.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        dxy(1)=-c(1);dxy(2)=-c(2);dxy(3)=-c(3);dxz(1)=b(1);dxz(2)=b(2);dxz(3)=b(3);
        for ii=1:3
            for jj=1:3,T11(ii,jj)=T11(ii,jj)+Weights2D(kt)*(rw(ii)*im(1,1)*rw(jj));
                       T12(ii,jj)=T12(ii,jj)+Weights2D(kt)*(rw(ii)*im(1,2)*dxy(jj)+rw(ii)*im(1,3)*dxz(jj));%
                       T13(ii,jj)=T13(ii,jj)+Weights2D(kt)*(rw(ii)*im(1,2)*wxy(jj)+rw(ii)*im(1,3)*wxz(jj));
                       T14(ii,jj)=T14(ii,jj)+Weights2D(kt)*(wxy(ii)*im(2,1)*rw(jj)+wxz(ii)*im(3,1)*rw(jj));
                       T15(ii,jj)=T15(ii,jj)+Weights2D(kt)*(wxy(ii)*im(2,2)*wxy(jj)+wxy(ii)*im(2,3)*wxz(jj)+wxz(ii)*im(3,2)*wxy(jj)+wxz(ii)*im(3,3)*wxz(jj));
                       T16(ii,jj)=T16(ii,jj)+Weights2D(kt)*(wxy(ii)*im(2,2)*dxy(jj)+wxy(ii)*im(2,3)*dxz(jj)+wxz(ii)*im(3,2)*dxy(jj)+wxz(ii)*im(3,3)*dxz(jj));
                       T17(ii,jj)=T17(ii,jj)+Weights2D(kt)*(dxy(ii)*im(2,1)*rw(jj)+dxz(ii)*im(3,1)*rw(jj));%
                       T18(ii,jj)=T18(ii,jj)+Weights2D(kt)*(dxy(ii)*im(2,2)*wxy(jj)+dxy(ii)*im(2,3)*wxz(jj)+dxz(ii)*im(3,2)*wxy(jj)+dxz(ii)*im(3,3)*wxz(jj));%
                       T19(ii,jj)=T19(ii,jj)+Weights2D(kt)*(dxy(ii)*im(2,2)*dxy(jj)+dxy(ii)*im(2,3)*dxz(jj)+dxz(ii)*im(3,2)*dxy(jj)+dxz(ii)*im(3,3)*dxz(jj));
                       T21(ii,jj)=T21(ii,jj)+Weights2D(kt)*(wwy(ii)*epsilon(2,2)*wwy(jj)+wwy(ii)*epsilon(2,3)*wwz(jj)+wwz(ii)*epsilon(3,2)*wwy(jj)+wwz(ii)*epsilon(3,3)*wwz(jj));
                       T22(ii,jj)=T22(ii,jj)+Weights2D(kt)*(wwy(ii)*epsilon(2,1)*simp(jj)+wwz(ii)*epsilon(3,1)*simp(jj));
                       T23(ii,jj)=T23(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(1,1)*simp(jj));
                       T24(ii,jj)=T24(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(1,2)*wwy(jj)+simp(ii)*epsilon(1,3)*wwz(jj));
             end
        end
    end
    for ii=1:3
        if(~isempty(edges(ii).OnLine) && edges(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
                if(toolboxModel.LineBoundaries(edges(ii).OnLine).Dispersive),cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param(FreqIndex);else,cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param;end,Tsge(ii,ii)=cond*edgeLength(ii);
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
        if(~isempty(vertices(ii).OnLine) && vertices(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
    end
        T11=T11*Ae;T12=T12*Ae;T13=T13*Ae;T14=T14*Ae;T15=T15*Ae;T16=T16*Ae;T17=T17*Ae;T18=T18*Ae;T19=T19*Ae;T21=T21*Ae;T22=T22*Ae;T23=T23*Ae;T24=T24*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-(omega^2)*si*sj*T21(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=si*sj*T11(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=1i*k0*waveImpedance*si*sj*Tsge(ii,jj);
                                   %------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=ej_E;L(counterL)=si*sj*(T14(ii,jj)-T13(ii,jj));
                                   %------------------------------------------------------------------------------------
                                   counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=-si*sj*T15(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=-(omega^2)*si*T22(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=si*T12(ii,jj);
                                   %-------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*T16(ii,jj);
            end
            if(vi_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=(omega^2)*T23(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=-T19(ii,jj);
            end
            if(vi_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=(omega^2)*sj*T24(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=-sj*T17(ii,jj);
                                   %-------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*T18(ii,jj);
            end
        end
    end
end
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_xBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==14),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};
                   medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;ksi=medium.Ksi;zi=medium.Zita;
                   frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==15),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};FrequencyIndex=varargin{15};medium=element.Medium2D;
                   if(medium.IsDispersive),epsilon=medium.Epsilon{FrequencyIndex};mu=medium.Mu{FrequencyIndex};ksi=medium.Ksi{FrequencyIndex};zi=medium.Zita{FrequencyIndex};
                   else,epsilon=medium.Epsilon;mu=medium.Mu;ksi=medium.Ksi;zi=medium.Zita;
                   end,frequency=toolboxModel.Frequency.Frequency(FrequencyIndex);
    end,im=(m0^-1)*mu^-1;epsilon=e0*epsilon;phi=ksi*im;psi=im*zi;th=ksi*im*zi;omega=2*pi*frequency;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];z=[vertices.Z];De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;b(2)=(z(3)-z(1))/De;b(3)=(z(1)-z(2))/De;c(1)=(y(3)-y(2))/De;c(2)=(y(1)-y(3))/De;c(3)=(y(2)-y(1))/De;   
    T11=zeros(3,3);T12=zeros(3,3);T13=zeros(3,3);T14=zeros(3,3);T15=zeros(3,3);T16=zeros(3,3);T17=zeros(3,3);T18=zeros(3,3);T19=zeros(3,3);
    T21=zeros(3,3);T22=zeros(3,3);T23=zeros(3,3);T24=zeros(3,3);T25=zeros(3,3);T26=zeros(3,3);
    T31=zeros(3,3);T32=zeros(3,3);T33=zeros(3,3);T34=zeros(3,3);T35=zeros(3,3);T36=zeros(3,3);
    T41=zeros(3,3);T42=zeros(3,3);T43=zeros(3,3);T44=zeros(3,3);T45=zeros(3,3);T46=zeros(3,3);T47=zeros(3,3);T48=zeros(3,3);
    Tsge=zeros(3,3);Tsgv=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwy(1)=simp(1)*b(2)-simp(2)*b(1);wwy(2)=simp(2)*b(3)-simp(3)*b(2);wwy(3)=simp(3)*b(1)-simp(1)*b(3);
        wwz(1)=simp(1)*c(2)-simp(2)*c(1);wwz(2)=simp(2)*c(3)-simp(3)*c(2);wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wxy(1)=-wwz(1);wxz(1)=wwy(1);wxy(2)=-wwz(2);wxz(2)=wwy(2);wxy(3)=-wwz(3);wxz(3)=wwy(3);
        wwy=wwy.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        dxy(1)=-c(1);dxy(2)=-c(2);dxy(3)=-c(3);dxz(1)=b(1);dxz(2)=b(2);dxz(3)=b(3);
        for ii=1:3
            for jj=1:3,T11(ii,jj)=T11(ii,jj)+Weights2D(kt)*(rw(ii)*im(1,1)*rw(jj));
                       T12(ii,jj)=T12(ii,jj)+Weights2D(kt)*(rw(ii)*im(1,2)*dxy(jj)+rw(ii)*im(1,3)*dxz(jj));%
                       T13(ii,jj)=T13(ii,jj)+Weights2D(kt)*(rw(ii)*im(1,2)*wxy(jj)+rw(ii)*im(1,3)*wxz(jj));
                       T14(ii,jj)=T14(ii,jj)+Weights2D(kt)*(wxy(ii)*im(2,1)*rw(jj)+wxz(ii)*im(3,1)*rw(jj));
                       T15(ii,jj)=T15(ii,jj)+Weights2D(kt)*(wxy(ii)*im(2,2)*wxy(jj)+wxy(ii)*im(2,3)*wxz(jj)+wxz(ii)*im(3,2)*wxy(jj)+wxz(ii)*im(3,3)*wxz(jj));
                       T16(ii,jj)=T16(ii,jj)+Weights2D(kt)*(wxy(ii)*im(2,2)*dxy(jj)+wxy(ii)*im(2,3)*dxz(jj)+wxz(ii)*im(3,2)*dxy(jj)+wxz(ii)*im(3,3)*dxz(jj));
                       T17(ii,jj)=T17(ii,jj)+Weights2D(kt)*(dxy(ii)*im(2,1)*rw(jj)+dxz(ii)*im(3,1)*rw(jj));%
                       T18(ii,jj)=T18(ii,jj)+Weights2D(kt)*(dxy(ii)*im(2,2)*wxy(jj)+dxy(ii)*im(2,3)*wxz(jj)+dxz(ii)*im(3,2)*wxy(jj)+dxz(ii)*im(3,3)*wxz(jj));%
                       T19(ii,jj)=T19(ii,jj)+Weights2D(kt)*(dxy(ii)*im(2,2)*dxy(jj)+dxy(ii)*im(2,3)*dxz(jj)+dxz(ii)*im(3,2)*dxy(jj)+dxz(ii)*im(3,3)*dxz(jj));
                       T21(ii,jj)=T21(ii,jj)+Weights2D(kt)*(rw(ii)*psi(1,2)*wwy(jj)+rw(ii)*psi(1,3)*wwz(jj));
                       T22(ii,jj)=T22(ii,jj)+Weights2D(kt)*(rw(ii)*psi(1,1)*simp(jj));
                       T23(ii,jj)=T23(ii,jj)+Weights2D(kt)*(dxy(ii)*psi(2,2)*wwy(jj)+dxy(ii)*psi(2,3)*wwz(jj)+dxz(ii)*psi(3,2)*wwy(jj)+dxz(ii)*psi(3,3)*wwz(jj));%
                       T24(ii,jj)=T24(ii,jj)+Weights2D(kt)*(dxy(ii)*psi(2,1)*simp(jj)+dxz(ii)*psi(3,1)*simp(jj));
                       T25(ii,jj)=T25(ii,jj)+Weights2D(kt)*(wxy(ii)*psi(2,2)*wwy(jj)+wxy(ii)*psi(2,3)*wwz(jj)+wxz(ii)*psi(3,2)*wwy(jj)+wxz(ii)*psi(3,3)*wwz(jj));
                       T26(ii,jj)=T26(ii,jj)+Weights2D(kt)*(wxy(ii)*psi(2,1)*simp(jj)+wxz(ii)*psi(3,1)*simp(jj));%
                       T31(ii,jj)=T31(ii,jj)+Weights2D(kt)*(wwy(ii)*phi(2,1)*rw(jj)+wwz(ii)*phi(3,1)*rw(jj));
                       T32(ii,jj)=T32(ii,jj)+Weights2D(kt)*(wwy(ii)*phi(2,2)*dxy(jj)+wwy(ii)*phi(2,3)*dxz(jj)+wwz(ii)*phi(3,2)*dxy(jj)+wwz(ii)*phi(3,3)*dxz(jj));
                       T33(ii,jj)=T33(ii,jj)+Weights2D(kt)*(simp(ii)*phi(1,1)*rw(jj));
                       T34(ii,jj)=T34(ii,jj)+Weights2D(kt)*(simp(ii)*phi(1,2)*dxy(jj)+simp(ii)*phi(1,3)*dxz(jj));
                       T35(ii,jj)=T35(ii,jj)+Weights2D(kt)*(wwy(ii)*phi(2,2)*wxy(jj)+wwy(ii)*phi(2,3)*wxz(jj)+wwz(ii)*phi(3,2)*wxy(jj)+wwz(ii)*phi(3,3)*wxz(jj));%
                       T36(ii,jj)=T36(ii,jj)+Weights2D(kt)*(simp(ii)*phi(1,2)*wxy(jj)+simp(ii)*phi(1,3)*wxz(jj));
                       T41(ii,jj)=T41(ii,jj)+Weights2D(kt)*(wwy(ii)*epsilon(2,2)*wwy(jj)+wwy(ii)*epsilon(2,3)*wwz(jj)+wwz(ii)*epsilon(3,2)*wwy(jj)+wwz(ii)*epsilon(3,3)*wwz(jj));
                       T42(ii,jj)=T42(ii,jj)+Weights2D(kt)*(wwy(ii)*th(2,2)*wwy(jj)+wwy(ii)*th(2,3)*wwz(jj)+wwz(ii)*th(3,2)*wwy(jj)+wwz(ii)*th(3,3)*wwz(jj));
                       T43(ii,jj)=T43(ii,jj)+Weights2D(kt)*(wwy(ii)*epsilon(2,1)*simp(jj)+wwz(ii)*epsilon(3,1)*simp(jj));
                       T44(ii,jj)=T44(ii,jj)+Weights2D(kt)*(wwy(ii)*th(2,1)*simp(jj)+wwz(ii)*th(3,1)*simp(jj));
                       T45(ii,jj)=T45(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(1,1)*simp(jj));
                       T46(ii,jj)=T46(ii,jj)+Weights2D(kt)*(simp(ii)*th(1,1)*simp(jj));
                       T47(ii,jj)=T47(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(1,2)*wwy(jj)+simp(ii)*epsilon(1,3)*wwz(jj));
                       T48(ii,jj)=T48(ii,jj)+Weights2D(kt)*(simp(ii)*th(1,2)*wwy(jj)+simp(ii)*th(1,3)*wwz(jj));
            end
        end
    end,n1=[1 2 3];n2=[2 3 1];
    for ii=1:3
        if(~isempty(edges(ii).OnLine) && edges(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
                if(toolboxModel.LineBoundaries(edges(ii).OnLine).Dispersive)
                    cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param(FreqIndex);
                else,cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param;
                end,Tsge(ii,ii)=cond*edgeLength(ii);
                    Tsgv(n1(ii),n1(ii))=cond*edgeLength(ii)/3;Tsgv(n2(ii),n2(ii))=cond*edgeLength(ii)/3;
                    Tsgv(n1(ii),n2(ii))=cond*edgeLength(ii)/6;Tsgv(n2(ii),n1(ii))=cond*edgeLength(ii)/3;
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
    end
    T11=T11*Ae;T12=T12*Ae;T13=T13*Ae;T14=T14*Ae;T15=T15*Ae;T16=T16*Ae;T17=T17*Ae;T18=T18*Ae;T19=T19*Ae;T21=T21*Ae;T22=T22*Ae;T23=T23*Ae;T24=T24*Ae;T25=T25*Ae;T26=T26*Ae;
    T31=T31*Ae;T32=T32*Ae;T33=T33*Ae;T34=T34*Ae;T35=T35*Ae;T36=T36*Ae;T41=T41*Ae;T42=T42*Ae;T43=T43*Ae;T44=T44*Ae;T45=T45*Ae;T46=T46*Ae;T47=T47*Ae;T48=T48*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=(omega^2)*si*sj*(T42(ii,jj)-T41(ii,jj));
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-1i*omega*si*sj*T31(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=1i*omega*si*sj*T21(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=si*sj*T11(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-1i*omega*si*sj*Tsge(ii,jj);
                                   %------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=ej_E;L(counterL)=si*sj*(T14(ii,jj)-T13(ii,jj)+1i*omega*T25(ii,jj)+1i*omega*T35(ii,jj));
                                   %------------------------------------------------------------------------------------
                                   counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=-si*sj*T15(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=(omega^2)*si*(T44(ii,jj)-T43(ii,jj));
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=-1i*omega*si*T32(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=1i*omega*si*T22(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=si*T12(ii,jj);
                                   %-------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*T16(ii,jj)+1i*omega*si*T26(ii,jj);
            end
            if(vi_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=(omega^2)*(T45(ii,jj)-T46(ii,jj));
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=1i*omega*T34(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=-1i*omega*T24(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=-T19(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=1i*omega*Tsgv(ii,jj);
            end
            if(vi_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=(omega^2)*sj*(T47(ii,jj)-T48(ii,jj));
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=1i*omega*sj*T33(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=-1i*omega*sj*T23(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=-sj*T17(ii,jj);
                                   %-------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*T18(ii,jj)-1i*omega*sj*T36(ii,jj);
            end
        end
    end
end
%---------------------------- y Boundaries --------------------------------
function [AssembledSystem] = VWE_Assembly_yBoundary(varargin)
    if(nargin==3),toolboxModel=varargin{1};AssembledSystem=varargin{2};boundaryIndices=varargin{3};
            N=AssembledSystem.DimEt+AssembledSystem.DimEn;
            IK=zeros(120*9*N,1);JK=zeros(120*9*N,1);K=zeros(120*9*N,1);IL=zeros(120*9*N,1);JL=zeros(120*9*N,1);L=zeros(120*9*N,1);
            IM=zeros(120*9*N,1);JM=zeros(120*9*N,1);M=zeros(120*9*N,1);counterK=0;counterL=0;counterM=0;
            for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
                for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                    switch medium.Type
                        case "Iso", [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_yBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,AssembledSystem,element);
                        case "Anis",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_yBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,AssembledSystem,element);
                        case "Bian",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_yBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,AssembledSystem,element);
                    end
                end
            end
            nonZerosK=nnz(IK);              nonZerosL=nnz(IL);               nonZerosM=nnz(IM);
            IK=IK(1:nonZerosK);             IL=IL(1:nonZerosL);              IM=IM(1:nonZerosM);
            JK=JK(1:nonZerosK);             JL=JL(1:nonZerosL);              JM=JM(1:nonZerosM);
            K=K(1:nonZerosK);               L=L(1:nonZerosL);                M=M(1:nonZerosM);
            K=sparse(IK,JK,K);              L=sparse(IL,JL,L);               M=sparse(IM,JM,M);

            Nv=AssembledSystem.DimEn;
            SA=sparse(2*Nv,2*Nv);           SB=sparse(2*Nv,2*Nv);
            SA(1:Nv,Nv+1:end)=-K;           SB(1:Nv,1:Nv)=-K;
            SA(1+Nv:end,1:Nv)=-K;           SB(Nv+1:end,Nv+1:end)=M;
            SA(1+Nv:end,1+Nv:end)=-L;

            AssembledSystem.Matrix_A=SA;    AssembledSystem.Matrix_B=SB;
    elseif(nargin==4),toolboxModel=varargin{1};AssembledSystem=varargin{2};boundaryIndices=varargin{3};FrequencyIndex=varargin{4};
        N=AssembledSystem.DimEt+AssembledSystem.DimEn;
            IK=zeros(120*9*N,1);JK=zeros(120*9*N,1);K=zeros(120*9*N,1);IL=zeros(120*9*N,1);JL=zeros(120*9*N,1);L=zeros(120*9*N,1);
            IM=zeros(120*9*N,1);JM=zeros(120*9*N,1);M=zeros(120*9*N,1);counterK=0;counterL=0;counterM=0;
            for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
                for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                    switch medium.Type
                        case "Iso", [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_yBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,AssembledSystem,element,FrequencyIndex);
                        case "Anis",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_yBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,AssembledSystem,element,FrequencyIndex);
                        case "Bian",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_yBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,AssembledSystem,element,FrequencyIndex);
                    end
                end
            end
            nonZerosK=nnz(IK);              nonZerosL=nnz(IL);               nonZerosM=nnz(IM);
            IK=IK(1:nonZerosK);             IL=IL(1:nonZerosL);              IM=IM(1:nonZerosM);
            JK=JK(1:nonZerosK);             JL=JL(1:nonZerosL);              JM=JM(1:nonZerosM);
            K=K(1:nonZerosK);               L=L(1:nonZerosL);                M=M(1:nonZerosM);
            K=sparse(IK,JK,K);              L=sparse(IL,JL,L);               M=sparse(IM,JM,M);

            Nv=AssembledSystem.DimEn;
            SA=sparse(2*Nv,2*Nv);           SB=sparse(2*Nv,2*Nv);
            SA(1:Nv,Nv+1:end)=-K;           SB(1:Nv,1:Nv)=-K;
            SA(1+Nv:end,1:Nv)=-K;           SB(Nv+1:end,Nv+1:end)=M;
            SA(1+Nv:end,1+Nv:end)=-L;

            AssembledSystem.Matrix_A{FrequencyIndex}=SA;    AssembledSystem.Matrix_B{FrequencyIndex}=SB;
    end
end
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_yBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==14),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};
                   medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;
                   frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==15),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};FrequencyIndex=varargin{15};medium=element.Medium2D;
                   if(medium.IsDispersive),epsilon=medium.Epsilon(FrequencyIndex);mu=medium.Mu(FrequencyIndex);
                   else,epsilon=medium.Epsilon;mu=medium.Mu;
                   end,frequency=toolboxModel.Frequency.Frequency(FrequencyIndex);
    end,im=(m0^-1)*mu^-1;epsilon=e0*epsilon;omega=2*pi*frequency;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];z=[vertices.Z];De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;b(2)=(z(3)-z(1))/De;b(3)=(z(1)-z(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;   
    T11=zeros(2,2);T15=zeros(2,2);T16=zeros(2,2);T18=zeros(2,2);T19=zeros(2,2);
    T21=zeros(2,2);T23=zeros(2,2);
    Tsge=zeros(2,2);Tsgv=zeros(2,2);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwz(1)=simp(1)*c(2)-simp(2)*c(1);wwz(2)=simp(2)*c(3)-simp(3)*c(2);wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wyx(1)=wwz(1);wyz(1)=-wwx(1);wyx(2)=wwz(2);wyz(2)=-wwx(2);wyx(3)=wwz(3);wyz(3)=-wwx(3);
        wwx=wwx.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        dyx(1)=c(1);dyx(2)=c(2);dyx(3)=c(3);dyz(1)=-b(1);dyz(2)=-b(2);dyz(3)=-b(3);
        for ii=1:3
            for jj=1:3,T11(ii,jj)=T11(ii,jj)+Weights2D(kt)*(rw(ii)*im*rw(jj));
                       T15(ii,jj)=T15(ii,jj)+Weights2D(kt)*(wyx(ii)*im*wyx(jj)+wyz(ii)*im*wyz(jj));
                       T16(ii,jj)=T16(ii,jj)+Weights2D(kt)*(wyx(ii)*im*dyx(jj)+wyz(ii)*im*dyz(jj));
                       T18(ii,jj)=T18(ii,jj)+Weights2D(kt)*(dyx(ii)*im*wyx(jj)+dyz(ii)*im*wyz(jj));
                       T19(ii,jj)=T19(ii,jj)+Weights2D(kt)*(dyx(ii)*im*dyx(jj)+dyz(ii)*im*dyz(jj));%
                       T21(ii,jj)=T21(ii,jj)+Weights2D(kt)*(wwx(ii)*epsilon*wwx(jj)+wwz(ii)*epsilon*wwz(jj));
                       T23(ii,jj)=T23(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon*simp(jj));
           end
        end
    end
    for ii=1:3
        if(~isempty(edges(ii).OnLine) && edges(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
                if(toolboxModel.LineBoundaries(edges(ii).OnLine).Dispersive),cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param(FreqIndex);else,cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param;end,Tsge(ii,ii)=cond*edgeLength(ii);
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
        if(~isempty(vertices(ii).OnLine) && vertices(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
    end
    T11=T11*Ae;T15=T15*Ae;T16=T16*Ae;T18=T18*Ae;T19=T19*Ae;T21=T21*Ae;T23=T23*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-(omega^2)*si*sj*T21(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=si*sj*T11(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=1i*k0*waveImpedance*si*sj*Tsge(ii,jj);
                                   %------------------------------------------------------------------------------------
                                   counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=-si*sj*T15(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*T16(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=(omega^2)*T23(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=-T19(ii,jj);
            end
            if(vi_E~=0 && ej_E~=0),counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*T18(ii,jj);end
        end
    end
end
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_yBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==14),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};
                   medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;
                   frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==15),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};FrequencyIndex=varargin{15};medium=element.Medium2D;
                   if(medium.IsDispersive),epsilon=medium.Epsilon{FrequencyIndex};mu=medium.Mu{FrequencyIndex};
                   else,epsilon=medium.Epsilon;mu=medium.Mu;
                   end,frequency=toolboxModel.Frequency.Frequency(FrequencyIndex);
    end,im=(m0^-1)*mu^-1;epsilon=e0*epsilon;omega=2*pi*frequency;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];z=[vertices.Z];De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;b(2)=(z(3)-z(1))/De;b(3)=(z(1)-z(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;   
    T11=zeros(2,2);T12=zeros(2,2);T13=zeros(2,2);T14=zeros(2,2);T15=zeros(2,2);T16=zeros(2,2);T17=zeros(2,2);T18=zeros(2,2);T19=zeros(2,2);
    T21=zeros(2,2);T22=zeros(2,2);T23=zeros(2,2);T24=zeros(2,2);
    Tsge=zeros(2,2);Tsgv=zeros(2,2);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwz(1)=simp(1)*c(2)-simp(2)*c(1);wwz(2)=simp(2)*c(3)-simp(3)*c(2);wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wyx(1)=wwz(1);wyz(1)=-wwx(1);wyx(2)=wwz(2);wyz(2)=-wwx(2);wyx(3)=wwz(3);wyz(3)=-wwx(3);
        wwx=wwx.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        dyx(1)=c(1);dyx(2)=c(2);dyx(3)=c(3);dyz(1)=-b(1);dyz(2)=-b(2);dyz(3)=-b(3);
        for ii=1:3
            for jj=1:3,T11(ii,jj)=T11(ii,jj)+Weights2D(kt)*(rw(ii)*im(2,2)*rw(jj));
                       T12(ii,jj)=T12(ii,jj)+Weights2D(kt)*(rw(ii)*im(2,3)*dyz(jj)+rw(ii)*im(2,1)*dyx(jj));
                       T13(ii,jj)=T13(ii,jj)+Weights2D(kt)*(rw(ii)*im(2,3)*wyz(jj)+rw(ii)*im(2,1)*wyx(jj));%
                       T14(ii,jj)=T14(ii,jj)+Weights2D(kt)*(wyx(ii)*im(1,2)*rw(jj)+wyz(ii)*im(3,2)*rw(jj));
                       T15(ii,jj)=T15(ii,jj)+Weights2D(kt)*(wyx(ii)*im(1,1)*wyx(jj)+wyx(ii)*im(1,3)*wyz(jj)+wyz(ii)*im(3,1)*wyx(jj)+wyz(ii)*im(3,3)*wyz(jj));
                       T16(ii,jj)=T16(ii,jj)+Weights2D(kt)*(wyx(ii)*im(1,1)*dyx(jj)+wyx(ii)*im(1,3)*dyz(jj)+wyz(ii)*im(3,1)*dyx(jj)+wyz(ii)*im(3,3)*dyz(jj));
                       T17(ii,jj)=T17(ii,jj)+Weights2D(kt)*(dyx(ii)*im(1,2)*rw(jj)+dyz(ii)*im(3,2)*rw(jj));%
                       T18(ii,jj)=T18(ii,jj)+Weights2D(kt)*(dyx(ii)*im(1,1)*wyx(jj)+dyx(ii)*im(1,3)*wyz(jj)+dyz(ii)*im(3,1)*wyx(jj)+dyz(ii)*im(3,3)*wyz(jj));
                       T19(ii,jj)=T19(ii,jj)+Weights2D(kt)*(dyx(ii)*im(1,1)*dyx(jj)+dyx(ii)*im(1,3)*dyz(jj)+dyz(ii)*im(3,1)*dyx(jj)+dyz(ii)*im(3,3)*dyz(jj));%
                       T21(ii,jj)=T21(ii,jj)+Weights2D(kt)*(wwx(ii)*epsilon(1,1)*wwx(jj)+wwx(ii)*epsilon(1,3)*wwz(jj)+wwz(ii)*epsilon(3,1)*wwx(jj)+wwz(ii)*epsilon(3,3)*wwz(jj));
                       T22(ii,jj)=T22(ii,jj)+Weights2D(kt)*(wwx(ii)*epsilon(1,2)*simp(jj)+wwz(ii)*epsilon(3,2)*simp(jj));
                       T23(ii,jj)=T23(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(2,2)*simp(jj));
                       T24(ii,jj)=T24(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(2,1)*wwx(jj)+simp(ii)*epsilon(2,3)*wwz(jj));                     
            end
        end
    end
    for ii=1:3
        if(~isempty(edges(ii).OnLine) && edges(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
                if(toolboxModel.LineBoundaries(edges(ii).OnLine).Dispersive),cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param(FreqIndex);else,cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param;end,Tsge(ii,ii)=cond*edgeLength(ii);
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
        if(~isempty(vertices(ii).OnLine) && vertices(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
    end
    T11=T11*Ae;T12=T12*Ae;T13=T13*Ae;T14=T14*Ae;T15=T15*Ae;T16=T16*Ae;T17=T17*Ae;T18=T18*Ae;T19=T19*Ae;T21=T21*Ae;T22=T22*Ae;T23=T23*Ae;T24=T24*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-(omega^2)*si*sj*T21(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=si*sj*T11(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=1i*k0*waveImpedance*si*sj*Tsge(ii,jj);
                                   %------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=ej_E;L(counterL)=si*sj*(T14(ii,jj)-T13(ii,jj));
                                   %------------------------------------------------------------------------------------
                                   counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=-si*sj*T15(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=-(omega^2)*si*T22(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=si*T12(ii,jj);
                                   %-------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*T16(ii,jj);
            end
            if(vi_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=(omega^2)*T23(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=-T19(ii,jj);
            end
            if(vi_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=(omega^2)*sj*T24(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=-sj*T17(ii,jj);
                                   %-------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*T18(ii,jj);
            end
        end
    end
end
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_yBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==14),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};
                   medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;ksi=medium.Ksi;zi=medium.Zita;
                   frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==15),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};FrequencyIndex=varargin{15};medium=element.Medium2D;
                   if(medium.IsDispersive),epsilon=medium.Epsilon{FrequencyIndex};mu=medium.Mu{FrequencyIndex};ksi=medium.Ksi{FrequencyIndex};zi=medium.Zita{FrequencyIndex};
                   else,epsilon=medium.Epsilon;mu=medium.Mu;ksi=medium.Ksi;zi=medium.Zita;
                   end,frequency=toolboxModel.Frequency.Frequency(FrequencyIndex);
    end,im=(m0^-1)*mu^-1;epsilon=e0*epsilon;phi=ksi*im;psi=im*zi;th=ksi*im*zi;omega=2*pi*frequency;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];z=[vertices.Z];De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;b(2)=(z(3)-z(1))/De;b(3)=(z(1)-z(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;   
    T11=zeros(2,2);T12=zeros(2,2);T13=zeros(2,2);T14=zeros(2,2);T15=zeros(2,2);T16=zeros(2,2);T17=zeros(2,2);T18=zeros(2,2);T19=zeros(2,2);
    T21=zeros(2,2);T22=zeros(2,2);T23=zeros(2,2);T24=zeros(2,2);T25=zeros(2,2);T26=zeros(2,2);
    T31=zeros(2,2);T32=zeros(2,2);T33=zeros(2,2);T34=zeros(2,2);T35=zeros(2,2);T36=zeros(2,2);
    T41=zeros(2,2);T42=zeros(2,2);T43=zeros(2,2);T44=zeros(2,2);T45=zeros(2,2);T46=zeros(2,2);T47=zeros(2,2);T48=zeros(2,2);
    Tsge=zeros(2,2);Tsgv=zeros(2,2);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwz(1)=simp(1)*c(2)-simp(2)*c(1);wwz(2)=simp(2)*c(3)-simp(3)*c(2);wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wyx(1)=wwz(1);wyz(1)=-wwx(1);wyx(2)=wwz(2);wyz(2)=-wwx(2);wyx(3)=wwz(3);wyz(3)=-wwx(3);
        wwx=wwx.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        dyx(1)=c(1);dyx(2)=c(2);dyx(3)=c(3);dyz(1)=-b(1);dyz(2)=-b(2);dyz(3)=-b(3);
        for ii=1:3
            for jj=1:3,T11(ii,jj)=T11(ii,jj)+Weights2D(kt)*(rw(ii)*im(2,2)*rw(jj));
                       T12(ii,jj)=T12(ii,jj)+Weights2D(kt)*(rw(ii)*im(2,3)*dyz(jj)+rw(ii)*im(2,1)*dyx(jj));
                       T13(ii,jj)=T13(ii,jj)+Weights2D(kt)*(rw(ii)*im(2,3)*wyz(jj)+rw(ii)*im(2,1)*wyx(jj));%
                       T14(ii,jj)=T14(ii,jj)+Weights2D(kt)*(wyx(ii)*im(1,2)*rw(jj)+wyz(ii)*im(3,2)*rw(jj));
                       T15(ii,jj)=T15(ii,jj)+Weights2D(kt)*(wyx(ii)*im(1,1)*wyx(jj)+wyx(ii)*im(1,3)*wyz(jj)+wyz(ii)*im(3,1)*wyx(jj)+wyz(ii)*im(3,3)*wyz(jj));
                       T16(ii,jj)=T16(ii,jj)+Weights2D(kt)*(wyx(ii)*im(1,1)*dyx(jj)+wyx(ii)*im(1,3)*dyz(jj)+wyz(ii)*im(3,1)*dyx(jj)+wyz(ii)*im(3,3)*dyz(jj));
                       T17(ii,jj)=T17(ii,jj)+Weights2D(kt)*(dyx(ii)*im(1,2)*rw(jj)+dyz(ii)*im(3,2)*rw(jj));%
                       T18(ii,jj)=T18(ii,jj)+Weights2D(kt)*(dyx(ii)*im(1,1)*wyx(jj)+dyx(ii)*im(1,3)*wyz(jj)+dyz(ii)*im(3,1)*wyx(jj)+dyz(ii)*im(3,3)*wyz(jj));
                       T19(ii,jj)=T19(ii,jj)+Weights2D(kt)*(dyx(ii)*im(1,1)*dyx(jj)+dyx(ii)*im(1,3)*dyz(jj)+dyz(ii)*im(3,1)*dyx(jj)+dyz(ii)*im(3,3)*dyz(jj));%
                       T21(ii,jj)=T21(ii,jj)+Weights2D(kt)*(rw(ii)*psi(2,1)*wwx(jj)+rw(ii)*psi(2,3)*wwz(jj));
                       T22(ii,jj)=T22(ii,jj)+Weights2D(kt)*(rw(ii)*psi(2,2)*simp(jj));
                       T23(ii,jj)=T23(ii,jj)+Weights2D(kt)*(dyx(ii)*psi(1,1)*wwx(jj)+dyx(ii)*psi(1,3)*wwz(jj)+dyz(ii)*psi(3,1)*wwx(jj)+dyz(ii)*psi(3,3)*wwz(jj));%
                       T24(ii,jj)=T24(ii,jj)+Weights2D(kt)*(dyx(ii)*psi(1,2)*simp(jj)+dyz(ii)*psi(3,2)*simp(jj));
                       T25(ii,jj)=T25(ii,jj)+Weights2D(kt)*(wyx(ii)*psi(1,1)*wwx(jj)+wyx(ii)*psi(1,3)*wwz(jj)+wyz(ii)*psi(3,1)*wwx(jj)+wyz(ii)*psi(3,3)*wwz(jj));
                       T26(ii,jj)=T26(ii,jj)+Weights2D(kt)*(wyx(ii)*psi(1,2)*simp(jj)+wyz(ii)*psi(3,2)*simp(jj));
                       T31(ii,jj)=T31(ii,jj)+Weights2D(kt)*(wwx(ii)*phi(1,2)*rw(jj)+wwz(ii)*phi(3,2)*rw(jj));
                       T32(ii,jj)=T32(ii,jj)+Weights2D(kt)*(wwx(ii)*phi(1,1)*dyx(jj)+wwx(ii)*phi(1,3)*dyz(jj)+wwz(ii)*phi(3,1)*dyx(jj)+wwz(ii)*phi(3,3)*dyz(jj));
                       T33(ii,jj)=T33(ii,jj)+Weights2D(kt)*(simp(ii)*phi(2,2)*rw(jj));
                       T34(ii,jj)=T34(ii,jj)+Weights2D(kt)*(simp(ii)*phi(2,1)*dyx(jj)+simp(ii)*phi(2,3)*dyz(jj));
                       T35(ii,jj)=T35(ii,jj)+Weights2D(kt)*(wwx(ii)*phi(1,1)*wyx(jj)+wwx(ii)*phi(1,3)*wyz(jj)+wwz(ii)*phi(3,1)*wyx(jj)+wwz(ii)*phi(3,3)*wyz(jj));
                       T36(ii,jj)=T36(ii,jj)+Weights2D(kt)*(simp(ii)*phi(2,1)*wyx(jj)+simp(ii)*phi(2,3)*wyz(jj));
                       T41(ii,jj)=T41(ii,jj)+Weights2D(kt)*(wwx(ii)*epsilon(1,1)*wwx(jj)+wwx(ii)*epsilon(1,3)*wwz(jj)+wwz(ii)*epsilon(3,1)*wwx(jj)+wwz(ii)*epsilon(3,3)*wwz(jj));
                       T42(ii,jj)=T42(ii,jj)+Weights2D(kt)*(wwx(ii)*th(1,1)*wwx(jj)+wwx(ii)*th(1,3)*wwz(jj)+wwz(ii)*th(3,1)*wwx(jj)+wwz(ii)*th(3,3)*wwz(jj));
                       T43(ii,jj)=T43(ii,jj)+Weights2D(kt)*(wwx(ii)*epsilon(1,2)*simp(jj)+wwz(ii)*epsilon(3,2)*simp(jj));
                       T44(ii,jj)=T44(ii,jj)+Weights2D(kt)*(wwx(ii)*th(1,2)*simp(jj)+wwz(ii)*th(3,2)*simp(jj));
                       T45(ii,jj)=T45(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(2,2)*simp(jj));
                       T46(ii,jj)=T46(ii,jj)+Weights2D(kt)*(simp(ii)*th(2,2)*simp(jj));
                       T47(ii,jj)=T47(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(2,1)*wwx(jj)+simp(ii)*epsilon(2,3)*wwz(jj));
                       T48(ii,jj)=T48(ii,jj)+Weights2D(kt)*(simp(ii)*th(2,1)*wwx(jj)+simp(ii)*th(2,3)*wwz(jj));
            end
        end


    end
    for ii=1:3
        if(~isempty(edges(ii).OnLine) && edges(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
                if(toolboxModel.LineBoundaries(edges(ii).OnLine).Dispersive),cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param(FreqIndex);else,cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param;end,Tsge(ii,ii)=cond*edgeLength(ii);
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
        if(~isempty(vertices(ii).OnLine) && vertices(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
    end
    T11=T11*Ae;T12=T12*Ae;T13=T13*Ae;T14=T14*Ae;T15=T15*Ae;T16=T16*Ae;T17=T17*Ae;T18=T18*Ae;T19=T19*Ae;T21=T21*Ae;T22=T22*Ae;T23=T23*Ae;T24=T24*Ae;T25=T25*Ae;T26=T26*Ae;
    T31=T31*Ae;T32=T32*Ae;T33=T33*Ae;T34=T34*Ae;T35=T35*Ae;T36=T36*Ae;T41=T41*Ae;T42=T42*Ae;T43=T43*Ae;T44=T44*Ae;T45=T45*Ae;T46=T46*Ae;T47=T47*Ae;T48=T48*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=(omega^2)*si*sj*(T42(ii,jj)-T41(ii,jj));
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-1i*omega*si*sj*T31(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=1i*omega*si*sj*T21(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=si*sj*T11(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=1i*k0*waveImpedance*si*sj*Tsge(ii,jj);
                                   %------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=ej_E;L(counterL)=si*sj*(T14(ii,jj)-T13(ii,jj)+1i*omega*T25(ii,jj)+1i*omega*T35(ii,jj));
                                   %------------------------------------------------------------------------------------
                                   counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=-si*sj*T15(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=(omega^2)*si*(T44(ii,jj)-T43(ii,jj));
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=-1i*omega*si*T32(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=1i*omega*si*T22(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=si*T12(ii,jj);
                                   %-------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*T16(ii,jj)+1i*omega*si*T26(ii,jj);
            end
            if(vi_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=(omega^2)*(T45(ii,jj)-T46(ii,jj));
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=1i*omega*T34(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=-1i*omega*T24(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=-T19(ii,jj);
            end
            if(vi_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=(omega^2)*sj*(T47(ii,jj)-T48(ii,jj));
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=1i*omega*sj*T33(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=-1i*omega*sj*T23(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=-sj*T17(ii,jj);
                                   %-------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*T18(ii,jj)-1i*omega*sj*T36(ii,jj);
            end
        end
    end
end
%---------------------------- z Boundaries --------------------------------
function [AssembledSystem] = VWE_Assembly_zBoundary(varargin)
    if(nargin==3),toolboxModel=varargin{1};AssembledSystem=varargin{2};boundaryIndices=varargin{3};
            N=AssembledSystem.DimEt+AssembledSystem.DimEn;
            IK=zeros(120*9*N,1);JK=zeros(120*9*N,1);K=zeros(120*9*N,1);IL=zeros(120*9*N,1);JL=zeros(120*9*N,1);L=zeros(120*9*N,1);
            IM=zeros(120*9*N,1);JM=zeros(120*9*N,1);M=zeros(120*9*N,1);counterK=0;counterL=0;counterM=0;
            for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
                for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                    switch medium.Type
                        case "Iso", [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_zBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,AssembledSystem,element);
                        case "Anis",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_zBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,AssembledSystem,element);
                        case "Bian",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_zBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,AssembledSystem,element);
                    end
                end
            end
            nonZerosK=nnz(IK);              nonZerosL=nnz(IL);               nonZerosM=nnz(IM);
            IK=IK(1:nonZerosK);             IL=IL(1:nonZerosL);              IM=IM(1:nonZerosM);
            JK=JK(1:nonZerosK);             JL=JL(1:nonZerosL);              JM=JM(1:nonZerosM);
            K=K(1:nonZerosK);               L=L(1:nonZerosL);                M=M(1:nonZerosM);
            K=sparse(IK,JK,K);              L=sparse(IL,JL,L);               M=sparse(IM,JM,M);

            Nv=AssembledSystem.DimEn;
            SA=sparse(2*Nv,2*Nv);           SB=sparse(2*Nv,2*Nv);
            SA(1:Nv,Nv+1:end)=-K;           SB(1:Nv,1:Nv)=-K;
            SA(1+Nv:end,1:Nv)=-K;           SB(Nv+1:end,Nv+1:end)=M;
            SA(1+Nv:end,1+Nv:end)=-L;

            AssembledSystem.Matrix_A=SA;    AssembledSystem.Matrix_B=SB;
    elseif(nargin==4),toolboxModel=varargin{1};AssembledSystem=varargin{2};boundaryIndices=varargin{3};FrequencyIndex=varargin{4};
        N=AssembledSystem.DimEt+AssembledSystem.DimEn;
            IK=zeros(120*9*N,1);JK=zeros(120*9*N,1);K=zeros(120*9*N,1);IL=zeros(120*9*N,1);JL=zeros(120*9*N,1);L=zeros(120*9*N,1);
            IM=zeros(120*9*N,1);JM=zeros(120*9*N,1);M=zeros(120*9*N,1);counterK=0;counterL=0;counterM=0;
            for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
                for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                    switch medium.Type
                        case "Iso", [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_zBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,AssembledSystem,element,FrequencyIndex);
                        case "Anis",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_zBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,AssembledSystem,element,FrequencyIndex);
                        case "Bian",[IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_zBoundary(IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM,toolboxModel,AssembledSystem,element,FrequencyIndex);
                    end
                end
            end
            nonZerosK=nnz(IK);              nonZerosL=nnz(IL);               nonZerosM=nnz(IM);
            IK=IK(1:nonZerosK);             IL=IL(1:nonZerosL);              IM=IM(1:nonZerosM);
            JK=JK(1:nonZerosK);             JL=JL(1:nonZerosL);              JM=JM(1:nonZerosM);
            K=K(1:nonZerosK);               L=L(1:nonZerosL);                M=M(1:nonZerosM);
            K=sparse(IK,JK,K);              L=sparse(IL,JL,L);               M=sparse(IM,JM,M);

            Nv=AssembledSystem.DimEn;
            SA=sparse(2*Nv,2*Nv);           SB=sparse(2*Nv,2*Nv);
            SA(1:Nv,Nv+1:end)=-K;           SB(1:Nv,1:Nv)=-K;
            SA(1+Nv:end,1:Nv)=-K;           SB(Nv+1:end,Nv+1:end)=M;
            SA(1+Nv:end,1+Nv:end)=-L;

            AssembledSystem.Matrix_A{FrequencyIndex}=SA;    AssembledSystem.Matrix_B{FrequencyIndex}=SB;
    end
end
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=IVWE_Assembly_zBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==14),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};
                   medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;
                   frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==15),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};FrequencyIndex=varargin{15};medium=element.Medium2D;
                   if(medium.IsDispersive),epsilon=medium.Epsilon(FrequencyIndex);mu=medium.Mu(FrequencyIndex);
                   else,epsilon=medium.Epsilon;mu=medium.Mu;
                   end,frequency=toolboxModel.Frequency.Frequency(FrequencyIndex);
    end,im=(m0^-1)*mu^-1;epsilon=e0*epsilon;omega=2*pi*frequency;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];x=[vertices.X];De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');Ae=abs(De)/2;
    b(1)=(y(2)-y(3))/De;b(2)=(y(3)-y(1))/De;b(3)=(y(1)-y(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;    
    T11=zeros(3,3);T15=zeros(3,3);T16=zeros(3,3);T18=zeros(3,3);T19=zeros(3,3);
    T21=zeros(3,3);T23=zeros(3,3);
    Tsge=zeros(3,3);Tsgv=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwy(1)=simp(1)*c(2)-simp(2)*c(1);wwy(2)=simp(2)*c(3)-simp(3)*c(2);wwy(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wzx(1)=-wwy(1);wzy(1)=wwx(1);wzx(2)=-wwy(2);wzy(2)=wwx(2);wzx(3)=-wwy(3);wzy(3)=wwx(3);
        wwx=wwx.*edgeLength;wwy=wwy.*edgeLength;rw=rw.*edgeLength;
        dzx(1)=-c(1);dzx(2)=-c(2);dzx(3)=-c(3);dzy(1)=b(1);dzy(2)=b(2);dzy(3)=b(3);
        for ii=1:3
            for jj=1:3,T11(ii,jj)=T11(ii,jj)+Weights2D(kt)*(rw(ii)*im*rw(jj));
                       T15(ii,jj)=T15(ii,jj)+Weights2D(kt)*(wzx(ii)*im*wzx(jj)+wzy(ii)*im*wzy(jj));
                       T16(ii,jj)=T16(ii,jj)+Weights2D(kt)*(wzx(ii)*im*dzx(jj)+wzy(ii)*im*dzy(jj));
                       T18(ii,jj)=T18(ii,jj)+Weights2D(kt)*(dzx(ii)*im*wzx(jj)+dzy(ii)*im*wzy(jj));
                       T19(ii,jj)=T19(ii,jj)+Weights2D(kt)*(dzx(ii)*im*dzx(jj)+dzy(ii)*im*dzy(jj));
                       T21(ii,jj)=T21(ii,jj)+Weights2D(kt)*(wwx(ii)*epsilon*wwx(jj)+wwy(ii)*epsilon*wwy(jj));
                       T23(ii,jj)=T23(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon*simp(jj));
            end
        end
    end
    for ii=1:3
        if(~isempty(edges(ii).OnLine) && edges(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
                if(toolboxModel.LineBoundaries(edges(ii).OnLine).Dispersive),cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param(FreqIndex);else,cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param;end,Tsge(ii,ii)=cond*edgeLength(ii);
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
        if(~isempty(vertices(ii).OnLine) && vertices(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
    end
    T11=T11*Ae;T15=T15*Ae;T16=T16*Ae;T18=T18*Ae;T19=T19*Ae;T21=T21*Ae;T23=T23*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-(omega^2)*si*sj*T21(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=si*sj*T11(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=1i*k0*waveImpedance*si*sj*Tsge(ii,jj);
                                   %------------------------------------------------------------------------------------
                                   counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=-si*sj*T15(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*T16(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=(omega^2)*T23(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=-T19(ii,jj);
            end
            if(vi_E~=0 && ej_E~=0),counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*T18(ii,jj);end
       end
    end
end
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=AVWE_Assembly_zBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==14),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};
                   medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;
                   frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==15),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};FrequencyIndex=varargin{15};medium=element.Medium2D;
                   if(medium.IsDispersive),epsilon=medium.Epsilon{FrequencyIndex};mu=medium.Mu{FrequencyIndex};
                   else,epsilon=medium.Epsilon;mu=medium.Mu;
                   end,frequency=toolboxModel.Frequency.Frequency(FrequencyIndex);
    end,im=(m0^-1)*mu^-1;epsilon=e0*epsilon;omega=2*pi*frequency;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];x=[vertices.X];De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');Ae=abs(De)/2;
    b(1)=(y(2)-y(3))/De;b(2)=(y(3)-y(1))/De;b(3)=(y(1)-y(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;     
    T11=zeros(3,3);T12=zeros(3,3);T13=zeros(3,3);T14=zeros(3,3);T15=zeros(3,3);T16=zeros(3,3);T17=zeros(3,3);T18=zeros(3,3);T19=zeros(3,3);
    T21=zeros(3,3);T22=zeros(3,3);T23=zeros(3,3);T24=zeros(3,3);
    Tsge=zeros(3,3);Tsgv=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwy(1)=simp(1)*c(2)-simp(2)*c(1);wwy(2)=simp(2)*c(3)-simp(3)*c(2);wwy(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wzx(1)=-wwy(1);wzy(1)=wwx(1);wzx(2)=-wwy(2);wzy(2)=wwx(2);wzx(3)=-wwy(3);wzy(3)=wwx(3);
        wwx=wwx.*edgeLength;wwy=wwy.*edgeLength;rw=rw.*edgeLength;
        dzx(1)=-c(1);dzx(2)=-c(2);dzx(3)=-c(3);dzy(1)=b(1);dzy(2)=b(2);dzy(3)=b(3);
        for ii=1:3
            for jj=1:3,T11(ii,jj)=T11(ii,jj)+Weights2D(kt)*(rw(ii)*im(3,3)*rw(jj));
                       T12(ii,jj)=T12(ii,jj)+Weights2D(kt)*(rw(ii)*im(3,2)*dzy(jj)+rw(ii)*im(3,1)*dzx(jj));
                       T13(ii,jj)=T13(ii,jj)+Weights2D(kt)*(rw(ii)*im(3,2)*wzy(jj)+rw(ii)*im(3,1)*wzx(jj));
                       T14(ii,jj)=T14(ii,jj)+Weights2D(kt)*(wzx(ii)*im(1,3)*rw(jj)+wzy(ii)*im(2,3)*rw(jj));
                       T15(ii,jj)=T15(ii,jj)+Weights2D(kt)*(wzx(ii)*im(1,1)*wzx(jj)+wzx(ii)*im(1,2)*wzy(jj)+wzy(ii)*im(2,1)*wzx(jj)+wzy(ii)*im(2,2)*wzy(jj));
                       T16(ii,jj)=T16(ii,jj)+Weights2D(kt)*(wzx(ii)*im(1,1)*dzx(jj)+wzx(ii)*im(1,2)*dzy(jj)+wzy(ii)*im(2,1)*dzx(jj)+wzy(ii)*im(2,2)*dzy(jj));
                       T17(ii,jj)=T17(ii,jj)+Weights2D(kt)*(dzx(ii)*im(1,3)*rw(jj)+dzy(ii)*im(2,3)*rw(jj));
                       T18(ii,jj)=T18(ii,jj)+Weights2D(kt)*(dzx(ii)*im(1,1)*wzx(jj)+dzx(ii)*im(1,2)*wzy(jj)+dzy(ii)*im(2,1)*wzx(jj)+dzy(ii)*im(2,2)*wzy(jj));
                       T19(ii,jj)=T19(ii,jj)+Weights2D(kt)*(dzx(ii)*im(1,1)*dzx(jj)+dzx(ii)*im(1,2)*dzy(jj)+dzy(ii)*im(2,1)*dzx(jj)+dzy(ii)*im(2,2)*dzy(jj));
                       T21(ii,jj)=T21(ii,jj)+Weights2D(kt)*(wwx(ii)*epsilon(1,1)*wwx(jj)+wwx(ii)*epsilon(1,2)*wwy(jj)+wwy(ii)*epsilon(2,1)*wwx(jj)+wwy(ii)*epsilon(2,2)*wwy(jj));
                       T22(ii,jj)=T22(ii,jj)+Weights2D(kt)*(wwx(ii)*epsilon(1,3)*simp(jj)+wwy(ii)*epsilon(2,3)*simp(jj));
                       T23(ii,jj)=T23(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(3,3)*simp(jj));
                       T24(ii,jj)=T24(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(3,1)*wwx(jj)+simp(ii)*epsilon(3,2)*wwy(jj));
            end
        end
    end
    for ii=1:3
        if(~isempty(edges(ii).OnLine) && edges(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
                if(toolboxModel.LineBoundaries(edges(ii).OnLine).Dispersive),cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param(FreqIndex);else,cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param;end,Tsge(ii,ii)=cond*edgeLength(ii);
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
        if(~isempty(vertices(ii).OnLine) && vertices(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
    end
    T11=T11*Ae;T12=T12*Ae;T13=T13*Ae;T14=T14*Ae;T15=T15*Ae;T16=T16*Ae;T17=T17*Ae;T18=T18*Ae;T19=T19*Ae;T21=T21*Ae;T22=T22*Ae;T23=T23*Ae;T24=T24*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-(omega^2)*si*sj*T21(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=si*sj*T11(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=1i*k0*waveImpedance*si*sj*Tsge(ii,jj);
                                   %------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=ej_E;L(counterL)=si*sj*(T14(ii,jj)-T13(ii,jj));
                                   %------------------------------------------------------------------------------------
                                   counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=-si*sj*T15(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=-(omega^2)*si*T22(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=si*T12(ii,jj);
                                   %-------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*T16(ii,jj);
            end
            if(vi_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=(omega^2)*T23(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=-T19(ii,jj);
            end
            if(vi_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=(omega^2)*sj*T24(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=-sj*T17(ii,jj);
                                   %-------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*T18(ii,jj);
            end
        end
    end
end
function [IK,JK,K,IL,JL,L,IM,JM,M,counterK,counterL,counterM]=BVWE_Assembly_zBoundary(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==14),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};
                   medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;ksi=medium.Ksi;zi=medium.Zita;
                   frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==15),IK=varargin{1};JK=varargin{2};K=varargin{3};IL=varargin{4};JL=varargin{5};L=varargin{6};IM=varargin{7};JM=varargin{8};M=varargin{9};counterK=varargin{10};counterL=varargin{11};counterM=varargin{10};
                   toolboxModel=varargin{13};element=varargin{14};FrequencyIndex=varargin{15};medium=element.Medium2D;
                   if(medium.IsDispersive),epsilon=medium.Epsilon{FrequencyIndex};mu=medium.Mu{FrequencyIndex};ksi=medium.Ksi{FrequencyIndex};zi=medium.Zita{FrequencyIndex};
                   else,epsilon=medium.Epsilon;mu=medium.Mu;ksi=medium.Ksi;zi=medium.Zita;
                   end,frequency=toolboxModel.Frequency.Frequency(FrequencyIndex);
    end,im=(m0^-1)*mu^-1;epsilon=e0*epsilon;phi=ksi*im;psi=im*zi;th=ksi*im*zi;omega=2*pi*frequency;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];x=[vertices.X];De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');Ae=abs(De)/2;
    b(1)=(y(2)-y(3))/De;b(2)=(y(3)-y(1))/De;b(3)=(y(1)-y(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;   
    T11=zeros(3,3);T12=zeros(3,3);T13=zeros(3,3);T14=zeros(3,3);T15=zeros(3,3);T16=zeros(3,3);T17=zeros(3,3);T18=zeros(3,3);T19=zeros(3,3);
    T21=zeros(3,3);T22=zeros(3,3);T23=zeros(3,3);T24=zeros(3,3);T25=zeros(3,3);T26=zeros(3,3);
    T31=zeros(3,3);T32=zeros(3,3);T33=zeros(3,3);T34=zeros(3,3);T35=zeros(3,3);T36=zeros(3,3);
    T41=zeros(3,3);T42=zeros(3,3);T43=zeros(3,3);T44=zeros(3,3);T45=zeros(3,3);T46=zeros(3,3);T47=zeros(3,3);T48=zeros(3,3);
    Tsge=zeros(3,3);Tsgv=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwy(1)=simp(1)*c(2)-simp(2)*c(1);wwy(2)=simp(2)*c(3)-simp(3)*c(2);wwy(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wzx(1)=-wwy(1);wzy(1)=wwx(1);wzx(2)=-wwy(2);wzy(2)=wwx(2);wzx(3)=-wwy(3);wzy(3)=wwx(3);
        wwx=wwx.*edgeLength;wwy=wwy.*edgeLength;rw=rw.*edgeLength;
        dzx(1)=-c(1);dzx(2)=-c(2);dzx(3)=-c(3);dzy(1)=b(1);dzy(2)=b(2);dzy(3)=b(3);
        for ii=1:3
            for jj=1:3,T11(ii,jj)=T11(ii,jj)+Weights2D(kt)*(rw(ii)*im(3,3)*rw(jj));
                       T12(ii,jj)=T12(ii,jj)+Weights2D(kt)*(rw(ii)*im(3,2)*dzy(jj)+rw(ii)*im(3,1)*dzx(jj));
                       T13(ii,jj)=T13(ii,jj)+Weights2D(kt)*(rw(ii)*im(3,2)*wzy(jj)+rw(ii)*im(3,1)*wzx(jj));
                       T14(ii,jj)=T14(ii,jj)+Weights2D(kt)*(wzx(ii)*im(1,3)*rw(jj)+wzy(ii)*im(2,3)*rw(jj));
                       T15(ii,jj)=T15(ii,jj)+Weights2D(kt)*(wzx(ii)*im(1,1)*wzx(jj)+wzx(ii)*im(1,2)*wzy(jj)+wzy(ii)*im(2,1)*wzx(jj)+wzy(ii)*im(2,2)*wzy(jj));
                       T16(ii,jj)=T16(ii,jj)+Weights2D(kt)*(wzx(ii)*im(1,1)*dzx(jj)+wzx(ii)*im(1,2)*dzy(jj)+wzy(ii)*im(2,1)*dzx(jj)+wzy(ii)*im(2,2)*dzy(jj));
                       T17(ii,jj)=T17(ii,jj)+Weights2D(kt)*(dzx(ii)*im(1,3)*rw(jj)+dzy(ii)*im(2,3)*rw(jj));
                       T18(ii,jj)=T18(ii,jj)+Weights2D(kt)*(dzx(ii)*im(1,1)*wzx(jj)+dzx(ii)*im(1,2)*wzy(jj)+dzy(ii)*im(2,1)*wzx(jj)+dzy(ii)*im(2,2)*wzy(jj));
                       T19(ii,jj)=T19(ii,jj)+Weights2D(kt)*(dzx(ii)*im(1,1)*dzx(jj)+dzx(ii)*im(1,2)*dzy(jj)+dzy(ii)*im(2,1)*dzx(jj)+dzy(ii)*im(2,2)*dzy(jj));
                       T21(ii,jj)=T21(ii,jj)+Weights2D(kt)*(rw(ii)*psi(3,1)*wwx(jj)+rw(ii)*psi(3,2)*wwy(jj));
                       T22(ii,jj)=T22(ii,jj)+Weights2D(kt)*(rw(ii)*psi(3,3)*simp(jj));
                       T23(ii,jj)=T23(ii,jj)+Weights2D(kt)*(dzx(ii)*psi(1,1)*wwx(jj)+dzx(ii)*psi(1,2)*wwy(jj)+dzy(ii)*psi(2,1)*wwx(jj)+dzy(ii)*psi(2,2)*wwy(jj));
                       T24(ii,jj)=T24(ii,jj)+Weights2D(kt)*(dzx(ii)*psi(1,3)*simp(jj)+dzy(ii)*psi(2,3)*simp(jj));
                       T25(ii,jj)=T25(ii,jj)+Weights2D(kt)*(wzx(ii)*psi(1,1)*wwx(jj)+wzx(ii)*psi(1,2)*wwy(jj)+wzy(ii)*psi(2,1)*wwx(jj)+wzy(ii)*psi(2,2)*wwy(jj));
                       T26(ii,jj)=T26(ii,jj)+Weights2D(kt)*(wzx(ii)*psi(1,3)*simp(jj)+wzy(ii)*psi(2,3)*simp(jj));
                       T31(ii,jj)=T31(ii,jj)+Weights2D(kt)*(wwx(ii)*phi(1,3)*rw(jj)+wwy(ii)*phi(2,3)*rw(jj));
                       T32(ii,jj)=T32(ii,jj)+Weights2D(kt)*(wwx(ii)*phi(1,1)*dzx(jj)+wwx(ii)*phi(1,2)*dzy(jj)+wwy(ii)*phi(2,1)*dzx(jj)+wwy(ii)*phi(2,2)*dzy(jj));
                       T33(ii,jj)=T33(ii,jj)+Weights2D(kt)*(simp(ii)*phi(3,3)*rw(jj));
                       T34(ii,jj)=T34(ii,jj)+Weights2D(kt)*(simp(ii)*phi(3,1)*dzx(jj)+simp(ii)*phi(3,2)*dzy(jj));
                       T35(ii,jj)=T35(ii,jj)+Weights2D(kt)*(wwx(ii)*phi(1,1)*wzx(jj)+wwx(ii)*phi(1,2)*wzy(jj)+wwy(ii)*phi(2,1)*wzx(jj)+wwy(ii)*phi(2,2)*wzy(jj));
                       T36(ii,jj)=T36(ii,jj)+Weights2D(kt)*(simp(ii)*phi(3,1)*wzx(jj)+simp(ii)*phi(3,2)*wzy(jj));
                       T41(ii,jj)=T41(ii,jj)+Weights2D(kt)*(wwx(ii)*epsilon(1,1)*wwx(jj)+wwx(ii)*epsilon(1,2)*wwy(jj)+wwy(ii)*epsilon(2,1)*wwx(jj)+wwy(ii)*epsilon(2,2)*wwy(jj));
                       T42(ii,jj)=T42(ii,jj)+Weights2D(kt)*(wwx(ii)*th(1,1)*wwx(jj)+wwx(ii)*th(1,2)*wwy(jj)+wwy(ii)*th(2,1)*wwx(jj)+wwy(ii)*th(2,2)*wwy(jj));
                       T43(ii,jj)=T43(ii,jj)+Weights2D(kt)*(wwx(ii)*epsilon(1,3)*simp(jj)+wwy(ii)*epsilon(2,3)*simp(jj));
                       T44(ii,jj)=T44(ii,jj)+Weights2D(kt)*(wwx(ii)*th(1,3)*simp(jj)+wwy(ii)*th(2,3)*simp(jj));
                       T45(ii,jj)=T45(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(3,3)*simp(jj));
                       T46(ii,jj)=T46(ii,jj)+Weights2D(kt)*(simp(ii)*th(3,3)*simp(jj));
                       T47(ii,jj)=T47(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(3,1)*wwx(jj)+simp(ii)*epsilon(3,2)*wwy(jj));
                       T48(ii,jj)=T48(ii,jj)+Weights2D(kt)*(simp(ii)*th(3,1)*wwx(jj)+simp(ii)*th(3,2)*wwy(jj));
            end
        end
    end
    for ii=1:3
        if(~isempty(edges(ii).OnLine) && edges(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
                if(toolboxModel.LineBoundaries(edges(ii).OnLine).Dispersive),cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param(FreqIndex);else,cond=toolboxModel.LineBoundaries(edges(ii).OnLine).Param;end,Tsge(ii,ii)=cond*edgeLength(ii);
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
        if(~isempty(vertices(ii).OnLine) && vertices(ii).OnLine~=0)
            if(toolboxModel.LineBoundaries(edges(ii).OnLine).Type=="GRA")
            elseif(toolboxModel.LineBoundaries(vertices(ii).OnLine).Type=="ABC")
            end
        end
    end
    T11=T11*Ae;T12=T12*Ae;T13=T13*Ae;T14=T14*Ae;T15=T15*Ae;T16=T16*Ae;T17=T17*Ae;T18=T18*Ae;T19=T19*Ae;T21=T21*Ae;T22=T22*Ae;T23=T23*Ae;T24=T24*Ae;T25=T25*Ae;T26=T26*Ae;
    T31=T31*Ae;T32=T32*Ae;T33=T33*Ae;T34=T34*Ae;T35=T35*Ae;T36=T36*Ae;T41=T41*Ae;T42=T42*Ae;T43=T43*Ae;T44=T44*Ae;T45=T45*Ae;T46=T46*Ae;T47=T47*Ae;T48=T48*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;vi_E=vertices(ii).IndexE;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;vj_E=vertices(jj).IndexE;
            if(ei_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=(omega^2)*si*sj*(T42(ii,jj)-T41(ii,jj));
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=-1i*omega*si*sj*T31(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=1i*omega*si*sj*T21(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=si*sj*T11(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=ej_E;K(counterK)=1i*k0*waveImpedance*si*sj*Tsge(ii,jj);
                                   %------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=ej_E;L(counterL)=si*sj*(T14(ii,jj)-T13(ii,jj)+1i*omega*T25(ii,jj)+1i*omega*T35(ii,jj));
                                   %------------------------------------------------------------------------------------
                                   counterM=counterM+1;IM(counterM)=ei_E;JM(counterM)=ej_E;M(counterM)=-si*sj*T15(ii,jj);
            end
            if(ei_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=(omega^2)*si*(T44(ii,jj)-T43(ii,jj));
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=-1i*omega*si*T32(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=1i*omega*si*T22(ii,jj);
                                   counterK=counterK+1;IK(counterK)=ei_E;JK(counterK)=vj_E;K(counterK)=si*T12(ii,jj);
                                   %-------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=ei_E;JL(counterL)=vj_E;L(counterL)=si*T16(ii,jj)+1i*omega*si*T26(ii,jj);
            end
            if(vi_E~=0 && vj_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=(omega^2)*(T45(ii,jj)-T46(ii,jj));
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=1i*omega*T34(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=-1i*omega*T24(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=vj_E;K(counterK)=-T19(ii,jj);
            end
            if(vi_E~=0 && ej_E~=0),counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=(omega^2)*sj*(T47(ii,jj)-T48(ii,jj));
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=1i*omega*sj*T33(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=-1i*omega*sj*T23(ii,jj);
                                   counterK=counterK+1;IK(counterK)=vi_E;JK(counterK)=ej_E;K(counterK)=-sj*T17(ii,jj);
                                   %-------------------------------------------------------------------------------------
                                   counterL=counterL+1;IL(counterL)=vi_E;JL(counterL)=ej_E;L(counterL)=sj*T18(ii,jj)-1i*omega*sj*T36(ii,jj);
            end
        end
    end
end