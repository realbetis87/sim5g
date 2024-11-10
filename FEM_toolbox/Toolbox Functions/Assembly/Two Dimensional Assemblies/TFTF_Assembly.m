%--------------------------------------------------------------------------
%{
              2D Transeverse Field Transverse Field Formulation 
                                (2D EH)

        Assembly for the 2D Transverse-Transverse-Field Method [1].
        Solution returns tangential components of intensity of electric
        field E and intensity of magnetic field H.(Used for non orthotropic media)

        [1] Zhu, Y. and Cangellaris, A.C. eds., 2006. Multigrid finite 
        element methods for electromagnetic field modeling (Vol. 28). 
        John Wiley & Sons. (9.4)
%}
%--------------------------------------------------------------------------
function [AssembledSystem] = TFTF_Assembly(toolboxModel,AssembledSystem,boundaryIndices),frequency=toolboxModel.Frequency;
    if(frequency.NF==1)
        if(isvector(boundaryIndices)),boundary=toolboxModel.Boundaries(boundaryIndices(1));
            switch boundary.Axis
                case  1,[AssembledSystem] = TFTF_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices,1);
                case -1,[AssembledSystem] = TFTF_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices,-1);
                case  2,[AssembledSystem] = TFTF_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices,1);
                case -2,[AssembledSystem] = TFTF_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices,-1);
                case  3,[AssembledSystem] = TFTF_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices,1);
                case -3,[AssembledSystem] = TFTF_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices,-1);
            end
        else,boundary=toolboxModel.Boundaries(boundaryIndices);
            switch boundary.Axis
                case  1,[AssembledSystem] = TFTF_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices,1);
                case -1,[AssembledSystem] = TFTF_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices,-1);
                case  2,[AssembledSystem] = TFTF_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices,1);
                case -2,[AssembledSystem] = TFTF_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices,-1);
                case  3,[AssembledSystem] = TFTF_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices,1);
                case -3,[AssembledSystem] = TFTF_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices,-1);
            end
        end
    else
        if(isvector(boundaryIndices)),boundary=toolboxModel.Boundaries(boundaryIndices(1));
            switch boundary.Axis
                case  1,for ii=1:frequency.NF,[AssembledSystem] = TFTF_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices,1,ii);end
                case -1,for ii=1:frequency.NF,[AssembledSystem] = TFTF_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices,-1,ii);end
                case  2,for ii=1:frequency.NF,[AssembledSystem] = TFTF_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices,1,ii);end
                case -2,for ii=1:frequency.NF,[AssembledSystem] = TFTF_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices,-1,ii);end
                case  3,for ii=1:frequency.NF,[AssembledSystem] = TFTF_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices,1,ii);end
                case -3,for ii=1:frequency.NF,[AssembledSystem] = TFTF_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices,-1,ii);end
            end
        else,boundary=toolboxModel.Boundaries(boundaryIndices);
            switch boundary.Axis
                case  1,for ii=1:frequency.NF,[AssembledSystem] = TFTF_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices,1,ii);end
                case -1,for ii=1:frequency.NF,[AssembledSystem] = TFTF_Assembly_xBoundary(toolboxModel,AssembledSystem,boundaryIndices,-1,ii);end
                case  2,for ii=1:frequency.NF,[AssembledSystem] = TFTF_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices,1,ii);end
                case -2,for ii=1:frequency.NF,[AssembledSystem] = TFTF_Assembly_yBoundary(toolboxModel,AssembledSystem,boundaryIndices,-1,ii);end
                case  3,for ii=1:frequency.NF,[AssembledSystem] = TFTF_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices,1,ii);end
                case -3,for ii=1:frequency.NF,[AssembledSystem] = TFTF_Assembly_zBoundary(toolboxModel,AssembledSystem,boundaryIndices,-1,ii);end
            end
        end
    end
end
%================== x Axis Propagation Assembly ===========================
function [AssembledSystem] = TFTF_Assembly_xBoundary(varargin)
    if(nargin==4),toolboxModel=varargin{1};AssembledSystem=varargin{2};boundaryIndices=varargin{3};nv=varargin{4};N=200;%N=AssembledSystem.DimEt+AssembledSystem.DimHt;
        II=zeros(120*9*N,1);JJ=zeros(120*9*N,1);SA=zeros(120*9*N,1);SB=zeros(120*9*N,1);counter=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case  "Iso",[II,JJ,SA,SB,counter] = ITFTF_Assembly_xBoundary(II,JJ,SA,SB,counter,toolboxModel,element,nv); 
                    case "Anis",[II,JJ,SA,SB,counter] = ATFTF_Assembly_xBoundary(II,JJ,SA,SB,counter,toolboxModel,element,nv); 
                end
            end
        end
    elseif(nargin==5),toolboxModel=varargin{1};AssembledSystem=varargin{2};boundaryIndices=varargin{3};nv=varargin{4};FrequencyIndex=varargin{5};N=AssembledSystem.DimEt+AssembledSystem.DimHt;
        II=zeros(120*9*N,1);JJ=zeros(120*9*N,1);SA=zeros(120*9*N,1);SB=zeros(120*9*N,1);counter=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case  "Iso",[II,JJ,SA,SB,counter] = ITFTF_Assembly_xBoundary(II,JJ,SA,SB,counter,toolboxModel,element,nv,FrequencyIndex);
                    case "Anis",[II,JJ,SA,SB,counter] = ATFTF_Assembly_xBoundary(II,JJ,SA,SB,counter,toolboxModel,element,nv,FrequencyIndex); 
                 end
            end
        end
    end
    nonZeros=nnz(II);II=II(1:nonZeros);JJ=JJ(1:nonZeros);SA=SA(1:nonZeros);SB=SB(1:nonZeros);
    SA=sparse(II,JJ,SA);SB=sparse(II,JJ,SB);AssembledSystem.Matrix_A=SA;AssembledSystem.Matrix_B=SB;
end
%--------------------------- Isotropic ------------------------------------
function [II,JJ,SA,SB,counter] = ITFTF_Assembly_xBoundary(varargin),ElectromagneticConstants;GaussianQuadratture2D;
    if(nargin==8),II=varargin{1};       SA=varargin{3};         counter=varargin{5};            element=varargin{7};
                  JJ=varargin{2};       SB=varargin{4};         toolboxModel=varargin{6};       nv=varargin{8};
        medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==9),II=varargin{1};   SA=varargin{3};         counter=varargin{5};            element=varargin{7};        FreqIndex=varargin{9};
                      JJ=varargin{2};   SB=varargin{4};         toolboxModel=varargin{6};       nv=varargin{8};      
        medium=element.Medium2D;
        if(medium.IsDispersive),epsilon=medium.Epsilon(FreqIndex);mu=medium.Mu(FreqIndex);frequency=toolboxModel.Frequency.Frequency(FreqIndex);
        else,epsilon=medium.Epsilon;mu=medium.Mu;frequency=toolboxModel.Frequency.Frequency;
        end
    end,omega=2*pi*frequency;k0=omega/c0;
    er=epsilon*eye(3,3);                                   mr=mu*eye(3,3);
    epsilon=er*e0;                                         mu=mr*m0;
    epsilonrc=[epsilon(2,2) epsilon(2,3);epsilon(3,2) epsilon(3,3);]-([epsilon(2,1);epsilon(3,1)]*[epsilon(1,2) epsilon(1,3)])/epsilon(1,1);   
    mrc=[mu(2,2) mu(2,3);mu(3,2) mu(3,3);]-([mu(2,1);mu(3,1)]*[mu(1,2) mu(1,3)])/mu(1,1);
    
    
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];z=[vertices.Z];De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;        c(1)=(y(3)-y(2))/De;
    b(2)=(z(3)-z(1))/De;        c(2)=(y(1)-y(3))/De;
    b(3)=(z(1)-z(2))/De;        c(3)=(y(2)-y(1))/De;
    %----------------------------------------------------------------------
    Se=zeros(3,3);Sm=zeros(3,3);Te=zeros(3,3);Tm=zeros(3,3);Pe=zeros(3,3);Pm=zeros(3,3);Qe=zeros(3,3);Qm=zeros(3,3);Hc=zeros(3,3);
    %----------------------------------------------------------------------
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wy(1)=simp(1)*b(2)-simp(2)*b(1);       wz(1)=simp(1)*c(2)-simp(2)*c(1);
        wy(2)=simp(2)*b(3)-simp(3)*b(2);       wz(2)=simp(2)*c(3)-simp(3)*c(2);
        wy(3)=simp(3)*b(1)-simp(1)*b(3);       wz(3)=simp(3)*c(1)-simp(1)*c(3);nv=-1;
        %------------------------------------------------------------------
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);          wxy(1)=-nv*wz(1); wxz(1)=nv*wy(1);       dxy(1)=-c(1);     dxz(1)=b(1);
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);          wxy(2)=-nv*wz(2); wxz(2)=nv*wy(2);       dxy(2)=-c(2);     dxz(2)=b(2);
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);          wxy(3)=-nv*wz(3); wxz(3)=nv*wy(3);       dxy(3)=-c(3);     dxz(3)=b(3);
        %------------------------------------------------------------------
        wy=wy.*edgeLength;wz=wz.*edgeLength;rw=rw.*edgeLength;
        wxy=wxy.*edgeLength;wxz=wxz.*edgeLength;
        for ii=1:3
            for jj=1:3
                Se(ii,jj) = Se(ii,jj) + Weights2D(kt) * (1/mu(1,1)) * (rw(ii) * rw(jj));
                Sm(ii,jj) = Sm(ii,jj) + Weights2D(kt) * (1/epsilon(1,1)) * (rw(ii) * rw(jj));
                Te(ii,jj) = Te(ii,jj) + Weights2D(kt) * (wy(ii)*epsilonrc(1,1)*wy(jj) +wy(ii)*epsilonrc(1,2)*wz(jj) ...
                                                        +wz(ii)*epsilonrc(2,1)*wy(jj) +wz(ii)*epsilonrc(2,2)*wz(jj));
                Tm(ii,jj) = Tm(ii,jj) + Weights2D(kt) * (wy(ii)*mrc(1,1)*wy(jj) + wy(ii)*mrc(1,2)*wz(jj) ...
                                                        +wz(ii)*mrc(2,1)*wy(jj) + wz(ii)*mrc(2,2)*wz(jj));
                Pe(ii,jj) = Pe(ii,jj) + Weights2D(kt) * (1/mu(1,1)) * (wy(ii) * mu(2,1) * rw(jj) + wz(ii) * mu(3,1) * rw(jj));
                Pm(ii,jj) = Pm(ii,jj) + Weights2D(kt) * (1/epsilon(1,1)) * (wy(ii) * epsilon(2,1) *rw(jj) + wz(ii) * epsilon(3,1) * rw(jj));
                Qe(ii,jj) = Qe(ii,jj) + Weights2D(kt) * (1/epsilon(1,1)) * (rw(ii) * epsilon(1,2) * wy(jj) + rw(ii) * epsilon(1,3) * wz(jj));
                Qm(ii,jj) = Qm(ii,jj) + Weights2D(kt) * (1/mu(1,1)) * (rw(ii) * mu(1,2) * wy(jj) + rw(ii) * mu(1,3) * wz(jj));

                Hc(ii,jj) = Hc(ii,jj) + Weights2D(kt) * (wy(ii) * wxy(jj) +wz(ii) * wxz(jj));
            end
        end
    end,Se=Se*Ae;Sm=Sm*Ae;Te=Te*Ae;Tm=Tm*Ae;Pe=Pe*Ae;Pm=Pm*Ae;Qe=Qe*Ae;Qm=Qm*Ae;Hc=Hc*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;
                                    II(counter)=ei_H;
                                    JJ(counter)=ej_H;
                                    SA(counter)=si*sj*Sm(ii,jj)-(omega^2)*si*sj*Tm(ii,jj);
            end
         
            if(ei_E~=0 && ej_E~=0),counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_E;
                                   SA(counter)=si*sj*Se(ii,jj)-(omega^2)*si*sj*Te(ii,jj);
            end            
            if(ei_H~=0 && ej_E~=0),counter=counter+1;
                                   II(counter)=ei_H;
                                   JJ(counter)=ej_E;
                                   SA(counter)=-1i*omega*si*sj*Pe(ii,jj)-1i*omega*si*sj*Qe(ii,jj);
                                   counter=counter+1;
                                   II(counter)=ei_H;
                                   JJ(counter)=ej_E;
                                   SB(counter)=omega*Hc(ii,jj);
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_H;
                                   SA(counter)=1i*omega*si*sj*Pm(ii,jj)+1i*omega*si*sj*Qm(ii,jj);
                                   counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_H;
                                   SB(counter)=-omega*Hc(ii,jj);
            end
        end
    end
end
%------------------ ------- Anisotropic -----------------------------------
function [II,JJ,SA,SB,counter] = ATFTF_Assembly_xBoundary(varargin),ElectromagneticConstants;GaussianQuadratture2D;
    if(nargin==8),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};element=varargin{7};nv=varargin{8};
        medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==9),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};element=varargin{7};nv=varargin{8};FreqIndex=varargin{9};
        medium=element.Medium2D;if(medium.IsDispersive),epsilon=medium.Epsilon{FreqIndex};mu=medium.Mu{FreqIndex};frequency=toolboxModel.Frequency.Frequency(FreqIndex);else,epsilon=medium.Epsilon;mu=medium.Mu;frequency=toolboxModel.Frequency.Frequency;end
    end
    mu=mu*m0;epsilon=epsilon*e0;omega=2*pi*frequency;
    mrc=[mu(2,2) mu(2,3);mu(3,2) mu(3,3);]-([mu(2,1);mu(3,1)]*[mu(1,2) mu(1,3)])/mu(1,1);
    epsilonrc=[epsilon(2,2) epsilon(2,3);epsilon(3,2) epsilon(3,3);]-([epsilon(2,1);epsilon(3,1)]*[epsilon(1,2) epsilon(1,3)])/epsilon(1,1);
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];z=[vertices.Z];De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;        c(1)=(y(3)-y(2))/De;
    b(2)=(z(3)-z(1))/De;        c(2)=(y(1)-y(3))/De;
    b(3)=(z(1)-z(2))/De;        c(3)=(y(2)-y(1))/De;
    %----------------------------------------------------------------------
    A1=zeros(3,3);A2=zeros(3,3);B1=zeros(3,3);B2=zeros(3,3);C11=zeros(3,3);C12=zeros(3,3);C21=zeros(3,3);C22=zeros(3,3);D=zeros(3,3);
    %----------------------------------------------------------------------
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wy(1)=simp(1)*b(2)-simp(2)*b(1);       wz(1)=simp(1)*c(2)-simp(2)*c(1);
        wy(2)=simp(2)*b(3)-simp(3)*b(2);       wz(2)=simp(2)*c(3)-simp(3)*c(2);
        wy(3)=simp(3)*b(1)-simp(1)*b(3);       wz(3)=simp(3)*c(1)-simp(1)*c(3);
        %------------------------------------------------------------------
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);          wxy(1)=-nv*wz(1); wxz(1)=nv*wy(1);       dxy(1)=-c(1);     dxz(1)=b(1);
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);          wxy(2)=-nv*wz(2); wxz(2)=nv*wy(2);       dxy(2)=-c(2);     dxz(2)=b(2);
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);          wxy(3)=-nv*wz(3); wxz(3)=nv*wy(3);       dxy(3)=-c(3);     dxz(3)=b(3);
        %------------------------------------------------------------------
        wy=wy.*edgeLength;wz=wz.*edgeLength;rw=rw.*edgeLength;
        for ii=1:3
            for jj=1:3
                A1(ii,jj) = A1(ii,jj) + Weights2D(kt) * (1/epsilon(1,1)) * ( rw(ii) * rw(jj) );
                A2(ii,jj) = A2(ii,jj) + Weights2D(kt) * ( wy(ii) * mrc(1,1) * wy(jj) + wy(ii) * mrc(1,2) * wz(jj) ...
                                                        + wz(ii) * mrc(2,1) * wy(jj) + wz(ii) * mrc(2,2) *wz(jj) );
                B1(ii,jj) = B1(ii,jj) + Weights2D(kt) * (1/mu(1,1)) * ( rw(ii) * rw(jj) );
                B2(ii,jj) = B2(ii,jj) + Weights2D(kt) * ( wy(ii) * epsilonrc(1,1) * wy(jj) + wy(ii) * epsilonrc(1,2) * wz(jj) ...
                                                        + wz(ii) * epsilonrc(2,1) * wy(jj) + wz(ii) * epsilonrc(2,2) * wz(jj) );
                C11(ii,jj) = C11(ii,jj) + Weights2D(kt) * (1/epsilon(1,1)) * (rw(ii) * epsilon(1,2) * wy(jj) + rw(ii) * epsilon(1,3) * wz(jj));
                C12(ii,jj) = C12(ii,jj) + Weights2D(kt) * (1/epsilon(1,1)) * (wy(ii) * epsilon(2,1) * rw(jj) + wz(ii) * epsilon(3,1) * rw(jj));
                C21(ii,jj) = C21(ii,jj) + Weights2D(kt) * (1/mu(1,1)) * (rw(ii) * mu(1,2) * wy(jj) + rw(ii) * mu(1,3) * wz(jj));
                C22(ii,jj) = C22(ii,jj) + Weights2D(kt) * (1/mu(1,1)) * (wy(ii) * mu(2,1) * rw(jj) + wz(ii) * mu(3,1) * rw(jj));
                D(ii,jj) = D(ii,jj) - Weights2D(kt) * ( wy(ii) * wxy(jj) + wz(ii) * wxz(jj) );
            end
        end
    end,A1 =A1*Ae;A2=A2*Ae;B1=B1*Ae;B2=B2*Ae;C11=C11*Ae;C12=C12*Ae;C21=C21*Ae;C22=C22*Ae;D=D*Ae;Dh=D';
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;
                                    II(counter)=ei_H;
                                    JJ(counter)=ej_H;
                                    SA(counter)=si*sj*A1(ii,jj)-(omega^2)*si*sj*A2(ii,jj);
            end
         
            if(ei_E~=0 && ej_E~=0),counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_E;
                                   SA(counter)=si*sj*B1(ii,jj)-(omega^2)*si*sj*B2(ii,jj);end
            
            if(ei_H~=0 && ej_E~=0),counter=counter+1;
                                   II(counter)=ei_H;
                                   JJ(counter)=ej_E;
                                   SA(counter)=-1i*omega*si*sj*C11(ii,jj)-1i*omega*si*sj*C12(ii,jj);
                                   counter=counter+1;
                                   II(counter)=ei_H;
                                   JJ(counter)=ej_E;
                                   SB(counter)=omega*D(ii,jj);
            
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_H;
                                   SA(counter)=1i*omega*si*sj*C21(ii,jj)+1i*omega*si*sj*C22(ii,jj);
                                   counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_H;
                                   SB(counter)=omega*Dh(ii,jj);
            end
        end
    end
end
%================== y Axis Propagation Assembly ===========================
function [AssembledSystem] = TFTF_Assembly_yBoundary(varargin)
    if(nargin==4),toolboxModel=varargin{1};AssembledSystem=varargin{2};boundaryIndices=varargin{3};nv=varargin{4};N=AssembledSystem.DimEt+AssembledSystem.DimHt;
        II=zeros(120*9*N,1);JJ=zeros(120*9*N,1);SA=zeros(120*9*N,1);SB=zeros(120*9*N,1);counter=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case  "Iso",[II,JJ,SA,SB,counter] = ITFTF_Assembly_yBoundary(II,JJ,SA,SB,counter,toolboxModel,element,nv); 
                    case "Anis",[II,JJ,SA,SB,counter] = ATFTF_Assembly_yBoundary(II,JJ,SA,SB,counter,toolboxModel,element,nv); 
                end
            end
        end
    elseif(nargin==5),toolboxModel=varargin{1};AssembledSystem=varargin{2};boundaryIndices=varargin{3};nv=varargin{4};FrequencyIndex=varargin{5};N=AssembledSystem.DimEt+AssembledSystem.DimHt;
        II=zeros(120*9*N,1);JJ=zeros(120*9*N,1);SA=zeros(120*9*N,1);SB=zeros(120*9*N,1);counter=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case  "Iso",[II,JJ,SA,SB,counter] = ITFTF_Assembly_yBoundary(II,JJ,SA,SB,counter,toolboxModel,element,nv,FrequencyIndex);
                    case "Anis",[II,JJ,SA,SB,counter] = ATFTF_Assembly_yBoundary(II,JJ,SA,SB,counter,toolboxModel,element,nv,FrequencyIndex); 
                 end
            end
        end
    end
    nonZeros=nnz(II);II=II(1:nonZeros);JJ=JJ(1:nonZeros);SA=SA(1:nonZeros);SB=SB(1:nonZeros);
    SA=sparse(II,JJ,SA);SB=sparse(II,JJ,SB);AssembledSystem.Matrix_A=SA;AssembledSystem.Matrix_B=SB;
end
%--------------------------- Isotropic ------------------------------------
function [II,JJ,SA,SB,counter] = ITFTF_Assembly_yBoundary(varargin),ElectromagneticConstants;GaussianQuadratture2D;
    if(nargin==8),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};element=varargin{7};nv=varargin{8};
        medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==9),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};element=varargin{7};nv=varargin{8};FreqIndex=varargin{9};
        medium=element.Medium2D;if(medium.IsDispersive),epsilon=medium.Epsilon(FreqIndex);mu=medium.Mu(FreqIndex);frequency=toolboxModel.Frequency.Frequency(FreqIndex);else,epsilon=medium.Epsilon;mu=medium.Mu;frequency=toolboxModel.Frequency.Frequency;end
    end
    er=epsilon*eye(3,3);mr=mu*epsilon*eye(3,3);epsilon=er;mu=mr;
    mu=mu*m0;epsilon=epsilon*e0;omega=2*pi*frequency;
    mrc=[mu(1,1) mu(1,3);mu(3,1) mu(3,3);]-([mu(1,2);mu(3,2)]*[mu(2,1) mu(2,3)])/mu(2,2);
    epsilonrc=[epsilon(1,1) epsilon(1,3);epsilon(3,1) epsilon(3,3);]-([epsilon(1,2);epsilon(3,2)]*[epsilon(2,1) epsilon(2,3)])/epsilon(2,2);
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];z=[vertices.Z];De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;        c(1)=(x(3)-x(2))/De;
    b(2)=(z(3)-z(1))/De;        c(2)=(x(1)-x(3))/De;
    b(3)=(z(1)-z(2))/De;        c(3)=(x(2)-x(1))/De;
    %----------------------------------------------------------------------
    A1=zeros(3,3);A2=zeros(3,3);B1=zeros(3,3);B2=zeros(3,3);C11=zeros(3,3);C12=zeros(3,3);C21=zeros(3,3);C22=zeros(3,3);D=zeros(3,3);
    %----------------------------------------------------------------------
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wx(1)=simp(1)*b(2)-simp(2)*b(1);       wz(1)=simp(1)*c(2)-simp(2)*c(1);
        wx(2)=simp(2)*b(3)-simp(3)*b(2);       wz(2)=simp(2)*c(3)-simp(3)*c(2);
        wx(3)=simp(3)*b(1)-simp(1)*b(3);       wz(3)=simp(3)*c(1)-simp(1)*c(3);
        %------------------------------------------------------------------
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);          wyx(1)=nv*wz(1); wyz(1)=-nv*wx(1);       dxy(1)=-c(1);     dxz(1)=b(1);
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);          wyx(2)=nv*wz(2); wyz(2)=-nv*wx(2);       dxy(2)=-c(2);     dxz(2)=b(2);
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);          wyx(3)=nv*wz(3); wyz(3)=-nv*wx(3);       dxy(3)=-c(3);     dxz(3)=b(3);
        %------------------------------------------------------------------
        wx=wx.*edgeLength;wz=wz.*edgeLength;rw=rw.*edgeLength;
        for ii=1:3
            for jj=1:3
                A1(ii,jj) = A1(ii,jj) + Weights2D(kt) * (1/epsilon(2,2)) * ( rw(ii) * rw(jj) );
                A2(ii,jj) = A2(ii,jj) + Weights2D(kt) * ( wx(ii) * mrc(1,1) * wx(jj) + wx(ii) * mrc(1,2) * wz(jj) ...
                                                        + wz(ii) * mrc(2,1) * wx(jj) + wz(ii) * mrc(2,2) *wz(jj) );
                B1(ii,jj) = B1(ii,jj) + Weights2D(kt) * (1/mu(2,2)) * ( rw(ii) * rw(jj) );
                B2(ii,jj) = B2(ii,jj) + Weights2D(kt) * ( wx(ii) * epsilonrc(1,1) * wx(jj) + wx(ii) * epsilonrc(1,2) * wz(jj) ...
                                                        + wz(ii) * epsilonrc(2,1) * wx(jj) + wz(ii) * epsilonrc(2,2) * wz(jj) );
                C11(ii,jj) = C11(ii,jj) + Weights2D(kt) * (1/epsilon(2,2)) * (rw(ii) * epsilon(2,1) * wx(jj) + rw(ii) * epsilon(2,3) * wz(jj));
                C12(ii,jj) = C12(ii,jj) + Weights2D(kt) * (1/epsilon(2,2)) * (wx(ii) * epsilon(1,2) * rw(jj) + wz(ii) * epsilon(3,2) * rw(jj));
                C21(ii,jj) = C21(ii,jj) + Weights2D(kt) * (1/mu(2,2)) * (rw(ii) * mu(2,1) * wx(jj) + rw(ii) * mu(2,3) * wz(jj));
                C22(ii,jj) = C22(ii,jj) + Weights2D(kt) * (1/mu(2,2)) * (wx(ii) * mu(1,2) * rw(jj) + wz(ii) * mu(3,2) * rw(jj));
                D(ii,jj) = D(ii,jj) - Weights2D(kt) * ( wx(ii) * wyx(jj) + wz(ii) * wyz(jj) );
            end
        end
    end,A1 =A1*Ae;A2=A2*Ae;B1=B1*Ae;B2=B2*Ae;C11=C11*Ae;C12=C12*Ae;C21=C21*Ae;C22=C22*Ae;D=D*Ae;Dh=D';
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;
                                    II(counter)=ei_H;
                                    JJ(counter)=ej_H;
                                    SA(counter)=si*sj*A1(ii,jj)-(omega^2)*si*sj*A2(ii,jj);
            end
         
            if(ei_E~=0 && ej_E~=0),counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_E;
                                   SA(counter)=si*sj*B1(ii,jj)-(omega^2)*si*sj*B2(ii,jj);end
            
            if(ei_H~=0 && ej_E~=0),counter=counter+1;
                                   II(counter)=ei_H;
                                   JJ(counter)=ej_E;
                                   SA(counter)=-1i*omega*si*sj*C11(ii,jj)-1i*omega*si*sj*C12(ii,jj);
                                   counter=counter+1;
                                   II(counter)=ei_H;
                                   JJ(counter)=ej_E;
                                   SB(counter)=omega*D(ii,jj);
            
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_H;
                                   SA(counter)=1i*omega*si*sj*C21(ii,jj)+1i*omega*si*sj*C22(ii,jj);
                                   counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_H;
                                   SB(counter)=omega*Dh(ii,jj);
            end
        end
    end
end
%------------------ ------- Anisotropic -----------------------------------
function [II,JJ,SA,SB,counter] = ATFTF_Assembly_yBoundary(varargin),ElectromagneticConstants;GaussianQuadratture2D;
    if(nargin==8),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};element=varargin{7};nv=varargin{8};
        medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==9),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};element=varargin{7};nv=varargin{8};FreqIndex=varargin{9};
        medium=element.Medium2D;if(medium.IsDispersive),epsilon=medium.Epsilon{FreqIndex};mu=medium.Mu{FreqIndex};frequency=toolboxModel.Frequency.Frequency(FreqIndex);else,epsilon=medium.Epsilon;mu=medium.Mu;frequency=toolboxModel.Frequency.Frequency;end
    end
    mu=mu*m0;epsilon=epsilon*e0;omega=2*pi*frequency;
    mrc=[mu(1,1) mu(1,3);mu(3,1) mu(3,3);]-([mu(1,2);mu(3,2)]*[mu(2,1) mu(2,3)])/mu(2,2);
    epsilonrc=[epsilon(1,1) epsilon(1,3);epsilon(3,1) epsilon(3,3);]-([epsilon(1,2);epsilon(3,2)]*[epsilon(2,1) epsilon(2,3)])/epsilon(2,2);
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];z=[vertices.Z];De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;        c(1)=(x(3)-x(2))/De;
    b(2)=(z(3)-z(1))/De;        c(2)=(x(1)-x(3))/De;
    b(3)=(z(1)-z(2))/De;        c(3)=(x(2)-x(1))/De;
    %----------------------------------------------------------------------
    A1=zeros(3,3);A2=zeros(3,3);B1=zeros(3,3);B2=zeros(3,3);C11=zeros(3,3);C12=zeros(3,3);C21=zeros(3,3);C22=zeros(3,3);D=zeros(3,3);
    %----------------------------------------------------------------------
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wx(1)=simp(1)*b(2)-simp(2)*b(1);       wz(1)=simp(1)*c(2)-simp(2)*c(1);
        wx(2)=simp(2)*b(3)-simp(3)*b(2);       wz(2)=simp(2)*c(3)-simp(3)*c(2);
        wx(3)=simp(3)*b(1)-simp(1)*b(3);       wz(3)=simp(3)*c(1)-simp(1)*c(3);
        %------------------------------------------------------------------
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);          wyx(1)=nv*wz(1); wyz(1)=-nv*wx(1);       dxy(1)=-c(1);     dxz(1)=b(1);
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);          wyx(2)=nv*wz(2); wyz(2)=-nv*wx(2);       dxy(2)=-c(2);     dxz(2)=b(2);
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);          wyx(3)=nv*wz(3); wyz(3)=-nv*wx(3);       dxy(3)=-c(3);     dxz(3)=b(3);
        %------------------------------------------------------------------
        wx=wx.*edgeLength;wz=wz.*edgeLength;rw=rw.*edgeLength;
        for ii=1:3
            for jj=1:3
                A1(ii,jj) = A1(ii,jj) + Weights2D(kt) * (1/epsilon(2,2)) * ( rw(ii) * rw(jj) );
                A2(ii,jj) = A2(ii,jj) + Weights2D(kt) * ( wx(ii) * mrc(1,1) * wx(jj) + wx(ii) * mrc(1,2) * wz(jj) ...
                                                        + wz(ii) * mrc(2,1) * wx(jj) + wz(ii) * mrc(2,2) *wz(jj) );
                B1(ii,jj) = B1(ii,jj) + Weights2D(kt) * (1/mu(2,2)) * ( rw(ii) * rw(jj) );
                B2(ii,jj) = B2(ii,jj) + Weights2D(kt) * ( wx(ii) * epsilonrc(1,1) * wx(jj) + wx(ii) * epsilonrc(1,2) * wz(jj) ...
                                                        + wz(ii) * epsilonrc(2,1) * wx(jj) + wz(ii) * epsilonrc(2,2) * wz(jj) );
                C11(ii,jj) = C11(ii,jj) + Weights2D(kt) * (1/epsilon(2,2)) * (rw(ii) * epsilon(2,1) * wx(jj) + rw(ii) * epsilon(2,3) * wz(jj));
                C12(ii,jj) = C12(ii,jj) + Weights2D(kt) * (1/epsilon(2,2)) * (wx(ii) * epsilon(1,2) * rw(jj) + wz(ii) * epsilon(3,2) * rw(jj));
                C21(ii,jj) = C21(ii,jj) + Weights2D(kt) * (1/mu(2,2)) * (rw(ii) * mu(2,1) * wx(jj) + rw(ii) * mu(2,3) * wz(jj));
                C22(ii,jj) = C22(ii,jj) + Weights2D(kt) * (1/mu(2,2)) * (wx(ii) * mu(1,2) * rw(jj) + wz(ii) * mu(3,2) * rw(jj));
                D(ii,jj) = D(ii,jj) - Weights2D(kt) * ( wx(ii) * wyx(jj) + wz(ii) * wyz(jj) );
            end
        end
    end,A1 =A1*Ae;A2=A2*Ae;B1=B1*Ae;B2=B2*Ae;C11=C11*Ae;C12=C12*Ae;C21=C21*Ae;C22=C22*Ae;D=D*Ae;Dh=D';
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;
                                    II(counter)=ei_H;
                                    JJ(counter)=ej_H;
                                    SA(counter)=si*sj*A1(ii,jj)-(omega^2)*si*sj*A2(ii,jj);
            end
         
            if(ei_E~=0 && ej_E~=0),counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_E;
                                   SA(counter)=si*sj*B1(ii,jj)-(omega^2)*si*sj*B2(ii,jj);end
            
            if(ei_H~=0 && ej_E~=0),counter=counter+1;
                                   II(counter)=ei_H;
                                   JJ(counter)=ej_E;
                                   SA(counter)=-1i*omega*si*sj*C11(ii,jj)-1i*omega*si*sj*C12(ii,jj);
                                   counter=counter+1;
                                   II(counter)=ei_H;
                                   JJ(counter)=ej_E;
                                   SB(counter)=omega*D(ii,jj);
            
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_H;
                                   SA(counter)=1i*omega*si*sj*C21(ii,jj)+1i*omega*si*sj*C22(ii,jj);
                                   counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_H;
                                   SB(counter)=omega*Dh(ii,jj);
            end
        end
    end
end
%================== z Axis Propagation Assembly ===========================
function [AssembledSystem] = TFTF_Assembly_zBoundary(varargin)
    if(nargin==4),toolboxModel=varargin{1};AssembledSystem=varargin{2};boundaryIndices=varargin{3};nv=varargin{4};N=AssembledSystem.DimEt+AssembledSystem.DimHt;
        II=zeros(120*9*N,1);JJ=zeros(120*9*N,1);SA=zeros(120*9*N,1);SB=zeros(120*9*N,1);counter=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case  "Iso",[II,JJ,SA,SB,counter] = ITFTF_Assembly_zBoundary(II,JJ,SA,SB,counter,toolboxModel,element,nv); 
                    case "Anis",[II,JJ,SA,SB,counter] = ATFTF_Assembly_zBoundary(II,JJ,SA,SB,counter,toolboxModel,element,nv); 
                end
            end
        end
    elseif(nargin==5),toolboxModel=varargin{1};AssembledSystem=varargin{2};boundaryIndices=varargin{3};nv=varargin{4};FrequencyIndex=varargin{5};N=AssembledSystem.DimEt+AssembledSystem.DimHt;
        II=zeros(120*9*N,1);JJ=zeros(120*9*N,1);SA=zeros(120*9*N,1);SB=zeros(120*9*N,1);counter=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case  "Iso",[II,JJ,SA,SB,counter] = ITFTF_Assembly_zBoundary(II,JJ,SA,SB,counter,toolboxModel,element,nv,FrequencyIndex);
                    case "Anis",[II,JJ,SA,SB,counter] = ATFTF_Assembly_zBoundary(II,JJ,SA,SB,counter,toolboxModel,element,nv,FrequencyIndex); 
                 end
            end
        end
    end
    nonZeros=nnz(II);II=II(1:nonZeros);JJ=JJ(1:nonZeros);SA=SA(1:nonZeros);SB=SB(1:nonZeros);
    SA=sparse(II,JJ,SA);SB=sparse(II,JJ,SB);AssembledSystem.Matrix_A=SA;AssembledSystem.Matrix_B=SB;
end
%--------------------------- Isotropic ------------------------------------
function [II,JJ,SA,SB,counter] = ITFTF_Assembly_zBoundary(varargin),ElectromagneticConstants;GaussianQuadratture2D;
    if(nargin==8),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};element=varargin{7};nv=varargin{8};
        medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==9),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};element=varargin{7};nv=varargin{8};FreqIndex=varargin{9};
        medium=element.Medium2D;if(medium.IsDispersive),epsilon=medium.Epsilon(FreqIndex);mu=medium.Mu(FreqIndex);frequency=toolboxModel.Frequency.Frequency(FreqIndex);else,epsilon=medium.Epsilon;mu=medium.Mu;frequency=toolboxModel.Frequency.Frequency;end
    end
    er=epsilon*eye(3,3);mr=mu*epsilon*eye(3,3);epsilon=er;mu=mr;
    mu=mu*m0;epsilon=epsilon*e0;omega=2*pi*frequency;
    mrc=[mu(1,1) mu(1,2);mu(2,1) mu(2,2);]-([mu(1,3);mu(2,3)]*[mu(3,1) mu(3,2)])/mu(3,3);
    epsilonrc=[epsilon(1,1) epsilon(1,2);epsilon(2,1) epsilon(2,2);]-([epsilon(1,3);epsilon(2,3)]*[epsilon(3,1) epsilon(3,2)])/epsilon(3,3);
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];y=[vertices.Y];De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');Ae=abs(De)/2;
    b(1)=(y(2)-y(3))/De;        c(1)=(x(3)-x(2))/De;
    b(2)=(y(3)-y(1))/De;        c(2)=(x(1)-x(3))/De;
    b(3)=(y(1)-y(2))/De;        c(3)=(x(2)-x(1))/De;
    %----------------------------------------------------------------------
    A1=zeros(3,3);A2=zeros(3,3);B1=zeros(3,3);B2=zeros(3,3);C11=zeros(3,3);C12=zeros(3,3);C21=zeros(3,3);C22=zeros(3,3);D=zeros(3,3);
    %----------------------------------------------------------------------
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wx(1)=simp(1)*b(2)-simp(2)*b(1);       wy(1)=simp(1)*c(2)-simp(2)*c(1);
        wx(2)=simp(2)*b(3)-simp(3)*b(2);       wy(2)=simp(2)*c(3)-simp(3)*c(2);
        wx(3)=simp(3)*b(1)-simp(1)*b(3);       wy(3)=simp(3)*c(1)-simp(1)*c(3);
        %------------------------------------------------------------------
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);          wzx(1)=-nv*wy(1); wzy(1)=nv*wx(1);       dxy(1)=-c(1);     dxz(1)=b(1);
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);          wzx(2)=-nv*wy(2); wzy(2)=nv*wx(2);       dxy(2)=-c(2);     dxz(2)=b(2);
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);          wzx(3)=-nv*wy(3); wzy(3)=nv*wx(3);       dxy(3)=-c(3);     dxz(3)=b(3);
        %------------------------------------------------------------------
        wx=wx.*edgeLength;wy=wy.*edgeLength;rw=rw.*edgeLength;
        for ii=1:3
            for jj=1:3
                A1(ii,jj) = A1(ii,jj) + Weights2D(kt) * (1/epsilon(3,3)) * ( rw(ii) * rw(jj) );
                A2(ii,jj) = A2(ii,jj) + Weights2D(kt) * ( wx(ii) * mrc(1,1) * wx(jj) + wx(ii) * mrc(1,2) * wy(jj) ...
                                                        + wy(ii) * mrc(2,1) * wx(jj) + wy(ii) * mrc(2,2) *wy(jj) );
                B1(ii,jj) = B1(ii,jj) + Weights2D(kt) * (1/mu(3,3)) * ( rw(ii) * rw(jj) );
                B2(ii,jj) = B2(ii,jj) + Weights2D(kt) * ( wx(ii) * epsilonrc(1,1) * wx(jj) + wx(ii) * epsilonrc(1,2) * wy(jj) ...
                                                        + wy(ii) * epsilonrc(2,1) * wx(jj) + wy(ii) * epsilonrc(2,2) * wy(jj) );
                C11(ii,jj) = C11(ii,jj) + Weights2D(kt) * (1/epsilon(3,3)) * (rw(ii) * epsilon(3,1) * wx(jj) + rw(ii) * epsilon(3,2) * wy(jj));
                C12(ii,jj) = C12(ii,jj) + Weights2D(kt) * (1/epsilon(3,3)) * (wx(ii) * epsilon(1,3) * rw(jj) + wy(ii) * epsilon(2,3) * rw(jj));
                C21(ii,jj) = C21(ii,jj) + Weights2D(kt) * (1/mu(3,3)) * (rw(ii) * mu(3,1) * wx(jj) + rw(ii) * mu(3,2) * wy(jj));
                C22(ii,jj) = C22(ii,jj) + Weights2D(kt) * (1/mu(3,3)) * (wx(ii) * mu(1,3) * rw(jj) + wy(ii) * mu(2,3) * rw(jj));
                D(ii,jj) = D(ii,jj) - Weights2D(kt) * ( wx(ii) * wzx(jj) + wy(ii) * wzy(jj) );
            end
        end
    end,A1 =A1*Ae;A2=A2*Ae;B1=B1*Ae;B2=B2*Ae;C11=C11*Ae;C12=C12*Ae;C21=C21*Ae;C22=C22*Ae;D=D*Ae;Dh=D';
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;
                                    II(counter)=ei_H;
                                    JJ(counter)=ej_H;
                                    SA(counter)=si*sj*A1(ii,jj)-(omega^2)*si*sj*A2(ii,jj);
            end
         
            if(ei_E~=0 && ej_E~=0),counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_E;
                                   SA(counter)=si*sj*B1(ii,jj)-(omega^2)*si*sj*B2(ii,jj);end
            
            if(ei_H~=0 && ej_E~=0),counter=counter+1;
                                   II(counter)=ei_H;
                                   JJ(counter)=ej_E;
                                   SA(counter)=-1i*omega*si*sj*C11(ii,jj)-1i*omega*si*sj*C12(ii,jj);
                                   counter=counter+1;
                                   II(counter)=ei_H;
                                   JJ(counter)=ej_E;
                                   SB(counter)=omega*D(ii,jj);
            
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_H;
                                   SA(counter)=1i*omega*si*sj*C21(ii,jj)+1i*omega*si*sj*C22(ii,jj);
                                   counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_H;
                                   SB(counter)=omega*Dh(ii,jj);
            end
        end
    end
end
%-------------------------- Anisotropic -----------------------------------
function [II,JJ,SA,SB,counter] = ATFTF_Assembly_zBoundary(varargin),ElectromagneticConstants;GaussianQuadratture2D;
    if(nargin==8),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};element=varargin{7};nv=varargin{8};
        medium=element.Medium2D;epsilon=medium.Epsilon;mu=medium.Mu;frequency=toolboxModel.Frequency.Frequency;
    elseif(nargin==9),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};element=varargin{7};nv=varargin{8};FreqIndex=varargin{9};
        medium=element.Medium2D;if(medium.IsDispersive),epsilon=medium.Epsilon{FreqIndex};mu=medium.Mu{FreqIndex};frequency=toolboxModel.Frequency.Frequency(FreqIndex);else,epsilon=medium.Epsilon;mu=medium.Mu;frequency=toolboxModel.Frequency.Frequency;end
    end
    mu=mu*m0;epsilon=epsilon*e0;omega=2*pi*frequency;
    mrc=[mu(1,1) mu(1,2);mu(2,1) mu(2,2);]-([mu(1,3);mu(2,3)]*[mu(3,1) mu(3,2)])/mu(3,3);
    epsilonrc=[epsilon(1,1) epsilon(1,2);epsilon(2,1) epsilon(2,2);]-([epsilon(1,3);epsilon(2,3)]*[epsilon(3,1) epsilon(3,2)])/epsilon(3,3);
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];y=[vertices.Y];De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');Ae=abs(De)/2;
    b(1)=(y(2)-y(3))/De;        c(1)=(x(3)-x(2))/De;
    b(2)=(y(3)-y(1))/De;        c(2)=(x(1)-x(3))/De;
    b(3)=(y(1)-y(2))/De;        c(3)=(x(2)-x(1))/De;
    %----------------------------------------------------------------------
    A1=zeros(3,3);A2=zeros(3,3);B1=zeros(3,3);B2=zeros(3,3);C11=zeros(3,3);C12=zeros(3,3);C21=zeros(3,3);C22=zeros(3,3);D=zeros(3,3);
    %----------------------------------------------------------------------
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wx(1)=simp(1)*b(2)-simp(2)*b(1);       wy(1)=simp(1)*c(2)-simp(2)*c(1);
        wx(2)=simp(2)*b(3)-simp(3)*b(2);       wy(2)=simp(2)*c(3)-simp(3)*c(2);
        wx(3)=simp(3)*b(1)-simp(1)*b(3);       wy(3)=simp(3)*c(1)-simp(1)*c(3);
        %------------------------------------------------------------------
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);          wzx(1)=-nv*wy(1); wzy(1)=nv*wx(1);       dxy(1)=-c(1);     dxz(1)=b(1);
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);          wzx(2)=-nv*wy(2); wzy(2)=nv*wx(2);       dxy(2)=-c(2);     dxz(2)=b(2);
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);          wzx(3)=-nv*wy(3); wzy(3)=nv*wx(3);       dxy(3)=-c(3);     dxz(3)=b(3);
        %------------------------------------------------------------------
        wx=wx.*edgeLength;wy=wy.*edgeLength;rw=rw.*edgeLength;
        for ii=1:3
            for jj=1:3
                A1(ii,jj) = A1(ii,jj) + Weights2D(kt) * (1/epsilon(3,3)) * ( rw(ii) * rw(jj) );
                A2(ii,jj) = A2(ii,jj) + Weights2D(kt) * ( wx(ii) * mrc(1,1) * wx(jj) + wx(ii) * mrc(1,2) * wy(jj) ...
                                                        + wy(ii) * mrc(2,1) * wx(jj) + wy(ii) * mrc(2,2) *wy(jj) );
                B1(ii,jj) = B1(ii,jj) + Weights2D(kt) * (1/mu(3,3)) * ( rw(ii) * rw(jj) );
                B2(ii,jj) = B2(ii,jj) + Weights2D(kt) * ( wx(ii) * epsilonrc(1,1) * wx(jj) + wx(ii) * epsilonrc(1,2) * wy(jj) ...
                                                        + wy(ii) * epsilonrc(2,1) * wx(jj) + wy(ii) * epsilonrc(2,2) * wy(jj) );
                C11(ii,jj) = C11(ii,jj) + Weights2D(kt) * (1/epsilon(3,3)) * (rw(ii) * epsilon(3,1) * wx(jj) + rw(ii) * epsilon(3,2) * wy(jj));
                C12(ii,jj) = C12(ii,jj) + Weights2D(kt) * (1/epsilon(3,3)) * (wx(ii) * epsilon(1,3) * rw(jj) + wy(ii) * epsilon(2,3) * rw(jj));
                C21(ii,jj) = C21(ii,jj) + Weights2D(kt) * (1/mu(3,3)) * (rw(ii) * mu(3,1) * wx(jj) + rw(ii) * mu(3,2) * wy(jj));
                C22(ii,jj) = C22(ii,jj) + Weights2D(kt) * (1/mu(3,3)) * (wx(ii) * mu(1,3) * rw(jj) + wy(ii) * mu(2,3) * rw(jj));
                D(ii,jj) = D(ii,jj) - Weights2D(kt) * ( wx(ii) * wzx(jj) + wy(ii) * wzy(jj) );
            end
        end
    end,A1 =A1*Ae;A2=A2*Ae;B1=B1*Ae;B2=B2*Ae;C11=C11*Ae;C12=C12*Ae;C21=C21*Ae;C22=C22*Ae;D=D*Ae;Dh=D';
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;
                                    II(counter)=ei_H;
                                    JJ(counter)=ej_H;
                                    SA(counter)=si*sj*A1(ii,jj)-(omega^2)*si*sj*A2(ii,jj);
            end
         
            if(ei_E~=0 && ej_E~=0),counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_E;
                                   SA(counter)=si*sj*B1(ii,jj)-(omega^2)*si*sj*B2(ii,jj);end
            
            if(ei_H~=0 && ej_E~=0),counter=counter+1;
                                   II(counter)=ei_H;
                                   JJ(counter)=ej_E;
                                   SA(counter)=-1i*omega*si*sj*C11(ii,jj)-1i*omega*si*sj*C12(ii,jj);
                                   counter=counter+1;
                                   II(counter)=ei_H;
                                   JJ(counter)=ej_E;
                                   SB(counter)=omega*D(ii,jj);
            
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_H;
                                   SA(counter)=1i*omega*si*sj*C21(ii,jj)+1i*omega*si*sj*C22(ii,jj);
                                   counter=counter+1;
                                   II(counter)=ei_E;
                                   JJ(counter)=ej_H;
                                   SB(counter)=omega*Dh(ii,jj);
            end
        end
    end
end