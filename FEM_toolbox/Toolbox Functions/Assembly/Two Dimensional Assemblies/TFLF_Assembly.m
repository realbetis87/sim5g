%--------------------------------------------------------------------------
%{
                 2D 1/2 Field Field Formulation
                           (2D 1/2 E-H)

    Assembly for the 2D 1/2 E-H Maxwelian Eigenmode solver for isotropic and
    orthotropic media. Solution returns tangential and normal to the domain 
    components of the intensity of the electric field E and intensity of the 
    magnetic field H.  
%}
%--------------------------------------------------------------------------
function [AssembledSystem] = TFLF_Assembly(toolboxModel,AssembledSystem,boundaryIndices),frequency=toolboxModel.Frequency;
     if(frequency.NF==1)
        if(isvector(boundaryIndices)),boundary=toolboxModel.Boundaries(boundaryIndices(1));freq=frequency.Frequency;
            switch abs(boundary.Axis)
                 case 1 
                    [MatA,MatB]=TFLF_Assembly_xBoundary(toolboxModel,boundaryIndices,freq);AssembledSystem.Matrix_A=MatA;AssembledSystem.Matrix_B=MatB;
                 case 2
                    [MatA,MatB]=TFLF_Assembly_yBoundary(toolboxModel,boundaryIndices,freq);AssembledSystem.Matrix_A=MatA;AssembledSystem.Matrix_B=MatB;
                 case 3
                    [MatA,MatB]=TFLF_Assembly_zBoundary(toolboxModel,boundaryIndices,freq);AssembledSystem.Matrix_A=MatA;AssembledSystem.Matrix_B=MatB;
             end
        else,boundary=toolboxModel.Boundaries(boundaryIndices);freq=frequency.Frequency;
            switch abs(boundary.Axis)
                 case 1 
                    [MatA,MatB]=TFLF_Assembly_xBoundary(toolboxModel,boundaryIndices,freq);AssembledSystem.Matrix_A=MatA;AssembledSystem.Matrix_B=MatB;
                 case 2
                    [MatA,MatB]=TFLF_Assembly_yBoundary(toolboxModel,boundaryIndices,freq);AssembledSystem.Matrix_A=MatA;AssembledSystem.Matrix_B=MatB;
                 case 3
                    [MatA,MatB]=TFLF_Assembly_zBoundary(toolboxModel,boundaryIndices,freq);AssembledSystem.Matrix_A=MatA;AssembledSystem.Matrix_B=MatB;
            end
        end
    else
       if(isvector(boundaryIndices)),boundary=toolboxModel.Boundaries(boundaryIndices(1));
            switch abs(boundary.Axis)
                case 1 
                    for ii=1:frequency.NF,freq=frequency.Frequency(ii);[MatA,MatB]=TFLF_Assembly_xBoundary(toolboxModel,boundaryIndices,freq,ii);AssembledSystem.Matrix_A{ii}=MatA;AssembledSystem.Matrix_B{ii}=MatB;end
                case 2
                    for ii=1:frequency.NF,freq=frequency.Frequency(ii);[MatA,MatB]=TFLF_Assembly_yBoundary(toolboxModel,boundaryIndices,freq,ii);AssembledSystem.Matrix_A{ii}=MatA;AssembledSystem.Matrix_B{ii}=MatB;end
                case 3
                    for ii=1:frequency.NF,freq=frequency.Frequency(ii);[MatA,MatB]=TFLF_Assembly_zBoundary(toolboxModel,boundaryIndices,freq,ii);AssembledSystem.Matrix_A{ii}=MatA;AssembledSystem.Matrix_B{ii}=MatB;end
             end
        else,boundary=toolboxModel.Boundaries(boundaryIndices);
            switch abs(boundary.Axis)
                case 1 
                    for ii=1:frequency.NF,freq=frequency.Frequency(ii);[MatA,MatB]=TFLF_Assembly_xBoundary(toolboxModel,boundaryIndices,freq,ii);AssembledSystem.Matrix_A{ii}=MatA;AssembledSystem.Matrix_B{ii}=MatB;end
                case 2
                    for ii=1:frequency.NF,freq=frequency.Frequency(ii);[MatA,MatB]=TFLF_Assembly_yBoundary(toolboxModel,boundaryIndices,freq,ii);AssembledSystem.Matrix_A{ii}=MatA;AssembledSystem.Matrix_B{ii}=MatB;end
                case 3
                    for ii=1:frequency.NF,freq=frequency.Frequency(ii);[MatA,MatB]=TFLF_Assembly_zBoundary(toolboxModel,boundaryIndices,freq,ii);AssembledSystem.Matrix_A{ii}=MatA;AssembledSystem.Matrix_B{ii}=MatB;end
            end
        end
    end
end
%---------------------- x Propagation Assembly ----------------------------
function [MatA,MatB] = TFLF_Assembly_xBoundary(varargin)
    if(nargin==3),toolboxModel=varargin{1};boundaryIndices=varargin{2};frequency=varargin{3};omega=2*pi*frequency;NumberOfElements=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));NumberOfElements=NumberOfElements+numel(boundary.Facets);end
        II = zeros(120*9*NumberOfElements,1); JJ = zeros(120*9*NumberOfElements,1);SA = zeros(120*9*NumberOfElements,1); SB =zeros(120*9*NumberOfElements,1);counter=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case "Iso"
                        [II,JJ,SA,SB,counter]=IsoAssembly_x(II,JJ,SA,SB,counter,toolboxModel,omega,element);
                    case "Anis"
                        [II,JJ,SA,SB,counter]=AnisAssembly_x(II,JJ,SA,SB,counter,toolboxModel,omega,element);
                end
            end
        end
    elseif(nargin==4),toolboxModel=varargin{1};boundaryIndices=varargin{2};frequency=varargin{3};FreqIndex=varargin{4};omega=2*pi*frequency;NumberOfElements=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));NumberOfElements=NumberOfElements+numel(boundary.Facets);end
        II = zeros(120*9*NumberOfElements,1); JJ = zeros(120*9*NumberOfElements,1);SA = zeros(120*9*NumberOfElements,1); SB =zeros(120*9*NumberOfElements,1);counter=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case "Iso"
                        [II,JJ,SA,SB,counter]=IsoAssembly_x(II,JJ,SA,SB,counter,toolboxModel,omega,element,FreqIndex);
                    case "Anis"
                        [II,JJ,SA,SB,counter]=AnisAssembly_x(II,JJ,SA,SB,counter,toolboxModel,omega,element,FreqIndex);
                end
            end
        end
    end
    NonZeros =nnz(II); SA=SA(1:NonZeros); II=II(1:NonZeros); JJ=JJ(1:NonZeros);SB=SB(1:NonZeros);
    MatA=sparse(II,JJ,SA);MatB=sparse(II,JJ,SB);
end
function [II,JJ,SA,SB,counter]=IsoAssembly_x(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==8),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};omega=varargin{7};element=varargin{8};
    elseif(nargin==9),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};omega=varargin{7};element=varargin{8};FreqIndex=varargin{9};
    end,medium=element.Medium2D;if(medium.IsDispersive),epsilon=medium.Epsilon(FreqIndex);mu=medium.Mu(FreqIndex);else,epsilon=medium.Epsilon;mu=medium.Mu;end,k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];z=[vertices.Z];De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;c(1)=(y(3)-y(2))/De;
    b(2)=(z(3)-z(1))/De;c(2)=(y(1)-y(3))/De;
    b(3)=(z(1)-z(2))/De;c(3)=(y(2)-y(1))/De;   
    Sh=zeros(3,3);Se=zeros(3,3);H=zeros(3,3);Th=zeros(3,3);Te=zeros(3,3);P=zeros(3,3);Tzh=zeros(3,3);Tze=zeros(3,3);Tse=zeros(3,3);Tsv=zeros(3,3);Tsh=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        %------------------------------------------------------------------
        wwy(1)=simp(1)*b(2)-simp(2)*b(1);
        wwy(2)=simp(2)*b(3)-simp(3)*b(2);
        wwy(3)=simp(3)*b(1)-simp(1)*b(3);
        wwz(1)=simp(1)*c(2)-simp(2)*c(1);
        wwz(2)=simp(2)*c(3)-simp(3)*c(2);
        wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);
        rw(2)=2*b(2)*c(3)-2*b(3)*c(2);
        rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wwy=wwy.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        %------------------------------------------------------------------
        wxy(1)=wwz(1);wxz(1)=-wwy(1);  dxy(1)=-c(1);dxz(1)=b(1);
        wxy(2)=wwz(2);wxz(2)=-wwy(2);  dxy(2)=-c(2);dxz(2)=b(2);
        wxy(3)=wwz(3);wxz(3)=-wwy(3);  dxy(3)=-c(3);dxz(3)=b(3);
        for ii=1:3
            for jj=1:3
                              Sh(ii,jj)=Sh(ii,jj)+Weights2D(kt)*(rw(ii)*(1/epsilon)*rw(jj));
                              Se(ii,jj)=Se(ii,jj)+Weights2D(kt)*(rw(ii)*(1/mu)*rw(jj));
                              H(ii,jj)=H(ii,jj)+Weights2D(kt)*(wwy(ii)*wxy(jj)+wwz(ii)*wxz(jj));
                              Th(ii,jj)=Th(ii,jj)+Weights2D(kt)*mu*(wwy(ii)*wwy(jj)+wwz(ii)*wwz(jj));
                              Te(ii,jj)=Te(ii,jj)+Weights2D(kt)*epsilon*(wwy(ii)*wwy(jj)+wwz(ii)*wwz(jj));
                              P(ii,jj)=P(ii,jj)+Weights2D(kt)*(simp(ii)*rw(jj));
                              Tze(ii,jj)=Tze(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon*simp(jj));
                              Tzh(ii,jj)=Tzh(ii,jj)+Weights2D(kt)*(simp(ii)*mu*simp(jj));
            end
        end
    end,Sh=Sh*Ae;Se=Se*Ae;H=H*Ae;Th=Th*Ae;Te=Te*Ae;P=P*Ae;Tzh=Tzh*Ae;Tze=Tze*Ae;%edgeSigns=[element.EdgeSigns];%edgeSigns=ones(6,1);
    nodes_1=[1 2 3];nodes_2=[2 3 1];
    for ii=1:3,edge=edges(ii);
        if(~isempty(edge.OnLine))
            if(edge.OnLine~=0),lineBoundary=toolboxModel.LineBoundaries(edge.OnLine);
                switch lineBoundary.Type
                    case "GRA"
                        if(lineBoundary.Dispersive),cond=lineBoundary.Param(freqIndex);else,cond=lineBoundary.Param;end
                        Tse(ii,ii)=-(k0^2)*mu*cond*element.EdgeSigns(ii)*edgeLength(ii)/(omega*e0);
                        Tsv(nodes_1(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/3;Tsv(nodes_2(ii),nodes_2(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/3;
                        Tsv(nodes_2(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/6;Tsv(nodes_2(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/6;
                    case "ABC"
                        Tse(ii,ii)=-1i*omega*m0*sqrt(epsilon/mu)*edgeLength(ii);
                        Tsh(ii,ii)=-1i*omega*e0*epsilon*sqrt(mu/epsilon)*edgeLength(ii);
                end
            end
        end
    end
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;vi_E=vertices(ii).IndexE;vi_H=vertices(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;vj_E=vertices(jj).IndexE;vj_H=vertices(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=-omega*m0*si*sj*Th(ii,jj);
                                   counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=(1/(omega*e0))*si*sj*Sh(ii,jj);
            end
            if(ei_E~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=omega*e0*si*sj*Te(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-(1/(omega*m0))*si*sj*Se(ii,jj);
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=si*sj*H(ii,jj);end
            if(ei_H~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_E;SB(counter)=si*sj*H(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=vj_E;SA(counter)=-1i*omega*e0*Tze(ii,jj);end
            if(vi_E~=0 && ej_H~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=ej_H;SA(counter)=sj*P(ii,jj);end
            if(vi_H~=0 && vj_H~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=vj_H;SA(counter)=1i*omega*m0*Tzh(ii,jj);end
            if(vi_H~=0 && ej_E~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=ej_E;SA(counter)=sj*P(ii,jj);end
        end
    end
    %{
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;vi_E=vertices(ii).IndexE;vi_H=vertices(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;vj_E=vertices(jj).IndexE;vj_H=vertices(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=-(k0^2)*si*sj*Th(ii,jj);%
                                   counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=si*sj*Sh(ii,jj)+si*sj*Tsh(ii,jj);
            end
            if(ei_E~=0 && ej_E~=0),%counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-(k0^2)*si*sj*Te(ii,jj)+si*sj*Tse(ii,jj);
                                   %counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=si*sj*Se(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=(k0^2)*si*sj*Te(ii,jj)+si*sj*Tse(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-si*sj*Se(ii,jj);
            end
            %if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=-(omega*m0)*si*sj*H(ii,jj);end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=(omega*m0)*si*sj*H(ii,jj);end
            if(ei_H~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_E;SB(counter)=(omega*e0)*si*sj*H(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=vj_E;SA(counter)=-1i*omega*e0*Tze(ii,jj)+Tsv(ii,jj);end
            if(vi_E~=0 && ej_H~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=ej_H;SA(counter)=sj*P(ii,jj);end
            if(vi_H~=0 && vj_H~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=vj_H;SA(counter)=1i*omega*m0*Tzh(ii,jj);end
            if(vi_H~=0 && ej_E~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=ej_E;SA(counter)=sj*P(ii,jj);end
        end
    end
    %}
end
function [II,JJ,SA,SB,counter]=AnisAssembly_x(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==8),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};omega=varargin{7};element=varargin{8};
    elseif(nargin==9),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};omega=varargin{7};element=varargin{8};FreqIndex=varargin{9};
    end,medium=element.Medium2D;if(medium.IsDispersive),epsilon=medium.Epsilon{FreqIndex};mu=medium.Mu{FreqIndex};else,epsilon=medium.Epsilon;mu=medium.Mu;end,k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];z=[vertices.Z];De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;b(2)=(z(3)-z(1))/De;b(3)=(z(1)-z(2))/De;c(1)=(y(3)-y(2))/De;c(2)=(y(1)-y(3))/De;c(3)=(y(2)-y(1))/De;   
    Sh=zeros(3,3);Se=zeros(3,3);H=zeros(3,3);Th=zeros(3,3);Te=zeros(3,3);P=zeros(3,3);Tzh=zeros(3,3);Tze=zeros(3,3);Tse=zeros(3,3);Tsv=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
         %------------------------------------------------------------------
        wwy(1)=simp(1)*b(2)-simp(2)*b(1);wwy(2)=simp(2)*b(3)-simp(3)*b(2);wwy(3)=simp(3)*b(1)-simp(1)*b(3);
        wwz(1)=simp(1)*c(2)-simp(2)*c(1);wwz(2)=simp(2)*c(3)-simp(3)*c(2);wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        %------------------------------------------------------------------
        wxy(1)=-wwz(1);wxz(1)=wwy(1);  dxy(1)=-c(1);dxz(1)=b(1);
        wxy(2)=-wwz(2);wxz(2)=wwy(2);  dxy(2)=-c(2);dxz(2)=b(2);
        wxy(3)=-wwz(3);wxz(3)=wwy(3);  dxy(3)=-c(3);dxz(3)=b(3);
        wwy=wwy.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        for ii=1:3,for jj=1:3,Sh(ii,jj)=Sh(ii,jj)+Weights2D(kt)*(rw(ii)*(1/epsilon(1,1))*rw(jj));
                              Se(ii,jj)=Se(ii,jj)+Weights2D(kt)*(rw(ii)*(1/mu(1,1))*rw(jj));
                              H(ii,jj)=H(ii,jj)+Weights2D(kt)*(wwy(ii)*wxy(jj)+wwz(ii)*wxz(jj));
                              Th(ii,jj)=Th(ii,jj)+Weights2D(kt)*(wwy(ii)*mu(2,2)*wwy(jj)+wwy(ii)*mu(2,3)*wwz(jj)+wwz(ii)*mu(3,3)*wwz(jj)+wwz(ii)*mu(3,2)*wwy(jj));
                              Te(ii,jj)=Te(ii,jj)+Weights2D(kt)*(wwy(ii)*epsilon(2,2)*wwy(jj)+wwy(ii)*epsilon(2,3)*wwz(jj)+wwz(ii)*epsilon(3,3)*wwz(jj)+wwz(ii)*epsilon(3,2)*wwy(jj));
                              P(ii,jj)=P(ii,jj)+Weights2D(kt)*(simp(ii)*rw(jj));
                              Tze(ii,jj)=Tze(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(1,1)*simp(jj));
                              Tzh(ii,jj)=Tzh(ii,jj)+Weights2D(kt)*(simp(ii)*mu(1,1)*simp(jj));
                   end
        end
    end,Sh=Sh*Ae;Se=Se*Ae;H=H*Ae;Th=Th*Ae;Te=Te*Ae;P=P*Ae;Tzh=Tzh*Ae;Tze=Tzh*Ae;
    nodes_1=[1 2 3];nodes_2=[2 3 1];
    for ii=1:3,edge=edges(ii);
        if(~isempty(edge.OnLine))
            if(edge.OnLine~=0),lineBoundary=toolboxModel.LineBoundaries(edge.OnLine);
                switch lineBoundary.Type
                    case "GRA"
                        if(lineBoundary.Dispersive),cond=lineBoundary.Param(freqIndex);else,cond=lineBoundary.Param;end
                        Tse(ii,ii)=-(k0^2)*mu*cond*element.EdgeSigns(ii)*edgeLength(ii)/(omega*e0);
                        Tsv(nodes_1(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/3;Tsv(nodes_2(ii),nodes_2(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/3;
                        Tsv(nodes_2(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/6;Tsv(nodes_2(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/6;
                    case "ABC"
                        Tse(ii,ii)=-1i*k0*sqrt(epsilon(2)*mu(2))*edgeLength(ii);
                end
            end
        end
    end
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;vi_E=vertices(ii).IndexE;vi_H=vertices(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;vj_E=vertices(jj).IndexE;vj_H=vertices(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=-(k0^2)*epsilon*si*sj*Th(ii,jj);%
                                   counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=si*sj*Sh(ii,jj);
            end
            if(ei_E~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-(k0^2)*mu*si*sj*Te(ii,jj)+si*sj*Tse(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=si*sj*Se(ii,jj);
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=-(omega*m0*mu)*si*sj*H(ii,jj);end
            if(ei_H~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_E;SB(counter)=(omega*e0*epsilon)*si*sj*H(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=vj_E;SA(counter)=-1i*omega*e0*Tze(ii,jj)+Tsv(ii,jj);end
            if(vi_E~=0 && ej_H~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=ej_H;SA(counter)=sj*P(ii,jj);end
            if(vi_H~=0 && vj_H~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=vj_H;SA(counter)=1i*omega*m0*Tzh(ii,jj);end
            if(vi_H~=0 && ej_E~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=ej_E;SA(counter)=sj*P(ii,jj);end
        end
    end
end
%---------------------- y Propagation Assembly ----------------------------
function [MatA,MatB] = TFLF_Assembly_yBoundary(varargin)
    if(nargin==3),toolboxModel=varargin{1};boundaryIndices=varargin{2};frequency=varargin{3};omega=2*pi*frequency;NumberOfElements=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));NumberOfElements=NumberOfElements+numel(boundary.Facets);end
        II = zeros(120*9*NumberOfElements,1); JJ = zeros(120*9*NumberOfElements,1);SA = zeros(120*9*NumberOfElements,1); SB =zeros(120*9*NumberOfElements,1);counter=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case "Iso"
                        [II,JJ,SA,SB,counter]=IsoAssembly_y(II,JJ,SA,SB,counter,toolboxModel,omega,element);
                    case "Anis"
                        [II,JJ,SA,SB,counter]=AnisAssembly_y(II,JJ,SA,SB,counter,toolboxModel,omega,element);
                end
            end
        end
    elseif(nargin==4),toolboxModel=varargin{1};boundaryIndices=varargin{2};frequency=varargin{3};FreqIndex=varargin{4};omega=2*pi*frequency;NumberOfElements=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));NumberOfElements=NumberOfElements+numel(boundary.Facets);end
        II = zeros(120*9*NumberOfElements,1); JJ = zeros(120*9*NumberOfElements,1);SA = zeros(120*9*NumberOfElements,1); SB =zeros(120*9*NumberOfElements,1);counter=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case "Iso"
                        [II,JJ,SA,SB,counter]=IsoAssembly_y(II,JJ,SA,SB,counter,toolboxModel,omega,element,FreqIndex);
                    case "Anis"
                        [II,JJ,SA,SB,counter]=AnisAssembly_y(II,JJ,SA,SB,counter,toolboxModel,omega,element,FreqIndex);
                end
            end
        end
    end
    NonZeros =nnz(II); SA=SA(1:NonZeros); II=II(1:NonZeros); JJ=JJ(1:NonZeros);SB=SB(1:NonZeros);
    MatA=sparse(II,JJ,SA);MatB=sparse(II,JJ,SB);
end
function [II,JJ,SA,SB,counter]=IsoAssembly_y(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==8),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};omega=varargin{7};element=varargin{8};
    elseif(nargin==9),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};omega=varargin{7};element=varargin{8};FreqIndex=varargin{9};
    end,medium=element.Medium2D;if(medium.IsDispersive),epsilon=medium.Epsilon(FreqIndex);mu=medium.Mu(FreqIndex);else,epsilon=medium.Epsilon;mu=medium.Mu;end,k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];z=[vertices.Z];De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;b(2)=(z(3)-z(1))/De;b(3)=(z(1)-z(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;   
    Sh=zeros(3,3);Se=zeros(3,3);H=zeros(3,3);Th=zeros(3,3);Te=zeros(3,3);P=zeros(3,3);Tzh=zeros(3,3);Tze=zeros(3,3);Tse=zeros(3,3);Tsv=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwz(1)=simp(1)*c(2)-simp(2)*c(1);wwz(2)=simp(2)*c(3)-simp(3)*c(2);wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wwy=wwy.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        wyx(1)=wwz(1);wyz(1)=-wwx(1);   dyx(1)=c(1);dyz(1)=-b(1);
        wyx(2)=wwz(2);wyz(2)=-wwx(2);   dyz(2)=c(2);dyz(2)=-b(2);
        wyx(3)=wwz(3);wyz(3)=-wwx(3);   dyz(3)=c(3);dyz(3)=-b(3);
        
        for ii=1:3,for jj=1:3,Sh(ii,jj)=Sh(ii,jj)+Weights2D(kt)*(rw(ii)*(1/epsilon)*rw(jj));
                              Se(ii,jj)=Se(ii,jj)+Weights2D(kt)*(rw(ii)*(1/mu)*rw(jj));
                              H(ii,jj)=H(ii,jj)+Weights2D(kt)*(wwx(ii)*wyx(jj)+wwz(ii)*wyz(jj));
                              Th(ii,jj)=Th(ii,jj)+Weights2D(kt)*mu*(wwx(ii)*wwx(jj)+wwz(ii)*wwz(jj));
                              Te(ii,jj)=Te(ii,jj)+Weights2D(kt)*epsilon*(wwx(ii)*wwx(jj)+wwz(ii)*wwz(jj));
                              P(ii,jj)=P(ii,jj)+Weights2D(kt)*(simp(ii)*rw(jj));
                              Tze(ii,jj)=Tze(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon*simp(jj));
                              Tzh(ii,jj)=Tzh(ii,jj)+Weights2D(kt)*(simp(ii)*mu*simp(jj));
                   end
        end
    end,Sh=Sh*Ae;Se=Se*Ae;H=H*Ae;Th=Th*Ae;Te=Te*Ae;P=P*Ae;Tzh=Tzh*Ae;Tze=Tze*Ae;%edgeSigns=[element.EdgeSigns];%edgeSigns=ones(6,1);
    nodes_1=[1 2 3];nodes_2=[2 3 1];
    for ii=1:3,edge=edges(ii);
        if(~isempty(edge.OnLine))
            if(edge.OnLine~=0),lineBoundary=toolboxModel.LineBoundaries(edge.OnLine);
                switch lineBoundary.Type
                    case "GRA"
                        if(lineBoundary.Dispersive),cond=lineBoundary.Param(freqIndex);else,cond=lineBoundary.Param;end
                        Tse(ii,ii)=-(k0^2)*mu*cond*element.EdgeSigns(ii)*edgeLength(ii)/(omega*e0);
                        Tsv(nodes_1(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/3;Tsv(nodes_2(ii),nodes_2(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/3;
                        Tsv(nodes_2(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/6;Tsv(nodes_2(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/6;
                    case "ABC"
                        Tse(ii,ii)=-1i*k0*sqrt(epsilon*mu)*edgeLength(ii);
                 end
            end
        end
    end
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;vi_E=vertices(ii).IndexE;vi_H=vertices(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;vj_E=vertices(jj).IndexE;vj_H=vertices(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=-(k0^2)*epsilon*si*sj*Th(ii,jj);%
                                   counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=si*sj*Sh(ii,jj);
            end
            if(ei_E~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-(k0^2)*mu*si*sj*Te(ii,jj)+si*sj*Tse(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=si*sj*Se(ii,jj);
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=-(omega*m0*mu)*si*sj*H(ii,jj);end
            if(ei_H~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_E;SB(counter)=(omega*e0*epsilon)*si*sj*H(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=vj_E;SA(counter)=-1i*omega*e0*Tze(ii,jj)+Tsv(ii,jj);end
            if(vi_E~=0 && ej_H~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=ej_H;SA(counter)=sj*P(ii,jj);end
            if(vi_H~=0 && vj_H~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=vj_H;SA(counter)=1i*omega*m0*Tzh(ii,jj);end
            if(vi_H~=0 && ej_E~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=ej_E;SA(counter)=sj*P(ii,jj);end
        end
    end
end
function [II,JJ,SA,SB,counter]=AnisAssembly_y(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==8),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};omega=varargin{7};element=varargin{8};
    elseif(nargin==9),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};omega=varargin{7};element=varargin{8};FreqIndex=varargin{9};
    end,medium=element.Medium2D;if(medium.IsDispersive),epsilon=medium.Epsilon{FreqIndex};mu=medium.Mu{FreqIndex};else,epsilon=medium.Epsilon;mu=medium.Mu;end,k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];z=[vertices.Z];De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;b(2)=(z(3)-z(1))/De;b(3)=(z(1)-z(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;   
    Sh=zeros(3,3);Se=zeros(3,3);H=zeros(3,3);Th=zeros(3,3);Te=zeros(3,3);P=zeros(3,3);Tzh=zeros(3,3);Tze=zeros(3,3);Tse=zeros(3,3);Tsv=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwz(1)=simp(1)*c(2)-simp(2)*c(1);wwz(2)=simp(2)*c(3)-simp(3)*c(2);wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wwx=wwx.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        wyx(1)=wwz(1);wyz(1)=-wwx(1);wyx(2)=wwz(2);wyz(2)=-wwx(2);wyx(3)=wwz(3);wyz(3)=-wwx(3);
        dyx(1)=c(1);dyz(1)=-b(1);dyz(2)=c(2);dyz(2)=-b(2);dyz(3)=c(3);dyz(3)=-b(3);
        for ii=1:3,for jj=1:3,Sh(ii,jj)=Sh(ii,jj)+Weights2D(kt)*(rw(ii)*(1/epsilon(2,2))*rw(jj));
                              Se(ii,jj)=Se(ii,jj)+Weights2D(kt)*(rw(ii)*(1/mu(2,2))*rw(jj));
                              H(ii,jj)=H(ii,jj)+Weights2D(kt)*(wwx(ii)*wyx(jj)+wwz(ii)*wyz(jj));
                              Th(ii,jj)=Th(ii,jj)+Weights2D(kt)*(wwx(ii)*mu(1,1)*wwx(jj)+wwx(ii)*mu(1,3)*wwz(jj)+wwz(ii)*mu(3,3)*wwz(jj)+wwz(ii)*mu(3,1)*wwx(jj));
                              Te(ii,jj)=Te(ii,jj)+Weights2D(kt)*(wwx(ii)*epsilon(1,1)*wwx(jj)+wwx(ii)*epsilon(1,3)*wwz(jj)+wwz(ii)*epsilon(3,3)*wwz(jj)+wwz(ii)*epsilon(3,1)*wwx(jj));
                              P(ii,jj)=P(ii,jj)+Weights2D(kt)*(simp(ii)*rw(jj));
                              Tze(ii,jj)=Tze(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon(2,2)*simp(jj));
                              Tzh(ii,jj)=Tzh(ii,jj)+Weights2D(kt)*(simp(ii)*mu(2,2)*simp(jj));
                   end
        end
    end,Sh=Sh*Ae;Se=Se*Ae;H=H*Ae;Th=Th*Ae;Te=Te*Ae;P=P*Ae;Tzh=Tzh*Ae;Tze=Tzh*Ae;
    nodes_1=[1 2 3];nodes_2=[2 3 1];
    for ii=1:3,edge=edges(ii);
        if(~isempty(edge.OnLine))
            if(edge.OnLine~=0),lineBoundary=toolboxModel.LineBoundaries(edge.OnLine);
                switch lineBoundary.Type
                    case "GRA"
                        if(lineBoundary.Dispersive),cond=lineBoundary.Param(freqIndex);else,cond=lineBoundary.Param;end
                        Tse(ii,ii)=-(k0^2)*mu*cond*element.EdgeSigns(ii)*edgeLength(ii)/(omega*e0);
                        Tsv(nodes_1(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/3;Tsv(nodes_2(ii),nodes_2(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/3;
                        Tsv(nodes_2(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/6;Tsv(nodes_2(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/6;
                    case "ABC"
                        Tse(ii,ii)=-1i*k0*sqrt(epsilon(2)*mu(2))*edgeLength(ii);
                end
            end
        end
    end
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;vi_E=vertices(ii).IndexE;vi_H=vertices(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;vj_E=vertices(jj).IndexE;vj_H=vertices(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=-(k0^2)*epsilon*si*sj*Th(ii,jj);%
                                   counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=si*sj*Sh(ii,jj);
            end
            if(ei_E~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-(k0^2)*mu*si*sj*Te(ii,jj)+si*sj*Tse(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=si*sj*Se(ii,jj);
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=-(omega*m0*mu)*si*sj*H(ii,jj);end
            if(ei_H~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_E;SB(counter)=(omega*e0*epsilon)*si*sj*H(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=vj_E;SA(counter)=-1i*omega*e0*Tze(ii,jj)+Tsv(ii,jj);end
            if(vi_E~=0 && ej_H~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=ej_H;SA(counter)=sj*P(ii,jj);end
            if(vi_H~=0 && vj_H~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=vj_H;SA(counter)=1i*omega*m0*Tzh(ii,jj);end
            if(vi_H~=0 && ej_E~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=ej_E;SA(counter)=sj*P(ii,jj);end
        end
    end
end
%---------------------- z Propagation Assembly ----------------------------
function [MatA,MatB] = TFLF_Assembly_zBoundary(varargin)
    if(nargin==3),toolboxModel=varargin{1};boundaryIndices=varargin{2};frequency=varargin{3};omega=2*pi*frequency;NumberOfElements=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));NumberOfElements=NumberOfElements+numel(boundary.Facets);end
        II = zeros(120*9*NumberOfElements,1); JJ = zeros(120*9*NumberOfElements,1);SA = zeros(120*9*NumberOfElements,1); SB =zeros(120*9*NumberOfElements,1);counter=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case "Iso"
                        [II,JJ,SA,SB,counter]=IsoAssembly_z(II,JJ,SA,SB,counter,toolboxModel,omega,element);
                    case "Anis"
                        [II,JJ,SA,SB,counter]=AnisAssembly_z(II,JJ,SA,SB,counter,toolboxModel,omega,element);
                end
            end
        end
    elseif(nargin==4),toolboxModel=varargin{1};boundaryIndices=varargin{2};frequency=varargin{3};omega=2*pi*frequency;FreqIndex=varargin{4};NumberOfElements=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));NumberOfElements=NumberOfElements+numel(boundary.Facets);end
        II = zeros(120*9*NumberOfElements,1); JJ = zeros(120*9*NumberOfElements,1);SA = zeros(120*9*NumberOfElements,1); SB =zeros(120*9*NumberOfElements,1);counter=0;
        for kk=1:numel(boundaryIndices),boundary=toolboxModel.Boundaries(boundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
                switch medium.Type
                    case "Iso"
                        [II,JJ,SA,SB,counter]=IsoAssembly_z(II,JJ,SA,SB,counter,toolboxModel,omega,element,FreqIndex);
                    case "Anis"
                        [II,JJ,SA,SB,counter]=AnisAssembly_z(II,JJ,SA,SB,counter,toolboxModel,omega,element,FreqIndex);
                end
            end
        end
    end
    NonZeros =nnz(II); SA=SA(1:NonZeros); II=II(1:NonZeros); JJ=JJ(1:NonZeros);SB=SB(1:NonZeros);
    MatA=sparse(II,JJ,SA);MatB=sparse(II,JJ,SB);
end
function [II,JJ,SA,SB,counter]=IsoAssembly_z(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==8),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};omega=varargin{7};element=varargin{8};
    elseif(nargin==9),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};omega=varargin{7};element=varargin{8};FreqIndex=varargin{9};
    end,medium=element.Medium2D;if(medium.IsDispersive),epsilon=medium.Epsilon(FreqIndex);mu=medium.Mu(FreqIndex);else,epsilon=medium.Epsilon;mu=medium.Mu;end,k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];y=[vertices.Y];De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');Ae=abs(De)/2;
    b(1)=(y(2)-y(3))/De;b(2)=(y(3)-y(1))/De;b(3)=(y(1)-y(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;   
    Sh=zeros(3,3);Se=zeros(3,3);H=zeros(3,3);Th=zeros(3,3);Te=zeros(3,3);P=zeros(3,3);Tzh=zeros(3,3);Tze=zeros(3,3);Tse=zeros(3,3);Tsv=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwy(1)=simp(1)*c(2)-simp(2)*c(1);wwy(2)=simp(2)*c(3)-simp(3)*c(2);wwy(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wwx=wwx.*edgeLength;wwy=wwy.*edgeLength;rw=rw.*edgeLength;
        wzx(1)=wwy(1);wzy(1)=-wwx(1);wzx(2)=wwy(2);wzy(2)=-wwx(2);wzx(3)=wwy(3);wzy(3)=-wwx(3);
        for ii=1:3,for jj=1:3,Sh(ii,jj)=Sh(ii,jj)+Weights2D(kt)*(rw(ii)*(1/epsilon)*rw(jj));
                              Se(ii,jj)=Se(ii,jj)+Weights2D(kt)*(rw(ii)*(1/mu)*rw(jj));
                              H(ii,jj)=H(ii,jj)+Weights2D(kt)*(wwx(ii)*wzx(jj)+wwy(ii)*wzy(jj));
                              Th(ii,jj)=Th(ii,jj)+Weights2D(kt)*mu*(wwx(ii)*wwx(jj)+wwy(ii)*wwy(jj));
                              Te(ii,jj)=Te(ii,jj)+Weights2D(kt)*epsilon*(wwx(ii)*wwx(jj)+wwy(ii)*wwy(jj));
                              P(ii,jj)=P(ii,jj)+Weights2D(kt)*(simp(ii)*rw(jj));
                              Tze(ii,jj)=Tze(ii,jj)+Weights2D(kt)*(simp(ii)*epsilon*simp(jj));
                              Tzh(ii,jj)=Tzh(ii,jj)+Weights2D(kt)*(simp(ii)*mu*simp(jj));
                   end
        end
    end,Sh=Sh*Ae;Se=Se*Ae;H=H*Ae;Th=Th*Ae;Te=Te*Ae;P=P*Ae;Tzh=Tzh*Ae;Tze=Tzh*Ae;
    nodes_1=[1 2 3];nodes_2=[2 3 1];
    for ii=1:3,edge=edges(ii);
        if(~isempty(edge.OnLine))
            if(edge.OnLine~=0),lineBoundary=toolboxModel.LineBoundaries(edge.OnLine);
                switch lineBoundary.Type
                    case "GRA"
                        if(lineBoundary.Dispersive),cond=lineBoundary.Param(freqIndex);else,cond=lineBoundary.Param;end
                        Tse(ii,ii)=-(k0^2)*mu*cond*element.EdgeSigns(ii)*edgeLength(ii)/(omega*e0);
                        Tsv(nodes_1(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/3;Tsv(nodes_2(ii),nodes_2(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/3;
                        Tsv(nodes_2(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/6;Tsv(nodes_2(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/6;
                    case "ABC"
                        Tse(ii,ii)=-1i*k0*sqrt(epsilon*mu)*edgeLength(ii);
                end
            end
        end
    end
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;vi_E=vertices(ii).IndexE;vi_H=vertices(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;vj_E=vertices(jj).IndexE;vj_H=vertices(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=-(k0^2)*epsilon*si*sj*Th(ii,jj);%
                                   counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=si*sj*Sh(ii,jj);
            end
            if(ei_E~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-(k0^2)*mu*si*sj*Te(ii,jj)+si*sj*Tse(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=si*sj*Se(ii,jj);
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=-(omega*m0*mu)*si*sj*H(ii,jj);end
            if(ei_H~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_E;SB(counter)=(omega*e0*epsilon)*si*sj*H(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=vj_E;SA(counter)=-1i*omega*e0*Tze(ii,jj)+Tsv(ii,jj);end
            if(vi_E~=0 && ej_H~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=ej_H;SA(counter)=sj*P(ii,jj);end
            if(vi_H~=0 && vj_H~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=vj_H;SA(counter)=1i*omega*m0*Tzh(ii,jj);end
            if(vi_H~=0 && ej_E~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=ej_E;SA(counter)=sj*P(ii,jj);end
        end
    end
end
function [II,JJ,SA,SB,counter]=AnisAssembly_z(varargin),GaussianQuadratture2D;ElectromagneticConstants;
    if(nargin==8),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};omega=varargin{7};element=varargin{8};
    elseif(nargin==9),II=varargin{1};JJ=varargin{2};SA=varargin{3};SB=varargin{4};counter=varargin{5};toolboxModel=varargin{6};omega=varargin{7};element=varargin{8};FreqIndex=varargin{9};
    end,medium=element.Medium2D;if(medium.IsDispersive),epsilon=medium.Epsilon{FreqIndex};mu=medium.Mu{FreqIndex};else,epsilon=medium.Epsilon;mu=medium.Mu;end,k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];y=[vertices.Y];De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');Ae=abs(De)/2;
    b(1)=(y(2)-y(3))/De;b(2)=(y(3)-y(1))/De;b(3)=(y(1)-y(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;   
    Sh=zeros(3,3);Se=zeros(3,3);H=zeros(3,3);Th=zeros(3,3);Te=zeros(3,3);P=zeros(3,3);Tzh=zeros(3,3);Tze=zeros(3,3);Tse=zeros(3,3);Tsv=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwy(1)=simp(1)*c(2)-simp(2)*c(1);wwy(2)=simp(2)*c(3)-simp(3)*c(2);wwy(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wwx=wwx.*edgeLength;wwy=wwy.*edgeLength;rw=rw.*edgeLength;
        wzx(1)=wwy(1);wzy(1)=-wwx(1);wzx(2)=wwy(2);wzy(2)=-wwx(2);wzx(3)=wwy(3);wzy(3)=-wwx(3);
        for ii=1:3,for jj=1:3,Sh(ii,jj)=Sh(ii,jj)+Weights2D(kt)*(rw(ii)*(1/epsilon(3,3))*rw(jj));
                              Se(ii,jj)=Se(ii,jj)+Weights2D(kt)*(rw(ii)*(1/mu(3,3))*rw(jj));
                              H(ii,jj)=H(ii,jj)+Weights2D(kt)*(wwx(ii)*wzx(jj)+wwy(ii)*wzy(jj));
                              Th(ii,jj)=Th(ii,jj)+Weights2D(kt)*(wwx(ii)*mu(1,1)*wwx(jj)+wwx(ii)*mu(1,2)*wwy(jj)+wwy(ii)*mu(2,2)*wwy(jj)+wwy(ii)*mu(2,1)*wwx(jj));
                              Te(ii,jj)=Te(ii,jj)+Weights2D(kt)*(wwx(ii)*epsilon(1,1)*wwx(jj)+wwx(ii)*epsilon(1,2)*wwy(jj)+wwy(ii)*epsilon(2,2)*wwy(jj)+wwy(ii)*epsilon(2,1)*wwx(jj));
                              P(ii,jj)=P(ii,jj)+Weights2D(kt)*(simp(ii)*rw(jj));
                              Tze(ii,jj)=Tze(ii,jj)+Weights2D(kt)*(simp(ii)*medium.Epsilon(3,3)*simp(jj));
                              Tzh(ii,jj)=Tzh(ii,jj)+Weights2D(kt)*(simp(ii)*mu(3,3)*simp(jj));
                   end
        end
    end,Sh=Sh*Ae;Se=Se*Ae;H=H*Ae;Th=Th*Ae;Te=Te*Ae;P=P*Ae;Tzh=Tzh*Ae;Tze=Tzh*Ae;
    nodes_1=[1 2 3];nodes_2=[2 3 1];
    for ii=1:3,edge=edges(ii);
        if(~isempty(edge.OnLine))
            if(edge.OnLine~=0),lineBoundary=toolboxModel.LineBoundaries(edge.OnLine);
                switch lineBoundary.Type
                    case "GRA"
                        if(lineBoundary.Dispersive),cond=lineBoundary.Param(freqIndex);else,cond=lineBoundary.Param;end
                        Tse(ii,ii)=-(k0^2)*mu*cond*element.EdgeSigns(ii)*edgeLength(ii)/(omega*e0);
                        Tsv(nodes_1(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/3;Tsv(nodes_2(ii),nodes_2(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/3;
                        Tsv(nodes_2(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/6;Tsv(nodes_2(ii),nodes_1(ii))=-1i*omega*e0*(-1i/(omega*e0))*cond*edgeLength(ii)/6;
                    case "ABC"
                        Tse(ii,ii)=-1i*k0*sqrt(epsilon(2)*mu(2))*edgeLength(ii);
                end
            end
        end
    end
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;vi_E=vertices(ii).IndexE;vi_H=vertices(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;vj_E=vertices(jj).IndexE;vj_H=vertices(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=-(k0^2)*epsilon*si*sj*Th(ii,jj);%
                                   counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=si*sj*Sh(ii,jj);
            end
            if(ei_E~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-(k0^2)*mu*si*sj*Te(ii,jj)+si*sj*Tse(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=si*sj*Se(ii,jj);
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=-(omega*m0*mu)*si*sj*H(ii,jj);end
            if(ei_H~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_E;SB(counter)=(omega*e0*epsilon)*si*sj*H(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=vj_E;SA(counter)=-1i*omega*e0*Tze(ii,jj)+Tsv(ii,jj);end
            if(vi_E~=0 && ej_H~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=ej_H;SA(counter)=sj*P(ii,jj);end
            if(vi_H~=0 && vj_H~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=vj_H;SA(counter)=1i*omega*m0*Tzh(ii,jj);end
            if(vi_H~=0 && ej_E~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=ej_E;SA(counter)=sj*P(ii,jj);end
        end
    end
end























%{
function [AssembledSystem] = TFLF_Assembly(toolboxModel,boundaryIndex),boundary=toolboxModel.Boundaries(boundaryIndex);frequency=toolboxModel.Frequency;
    AssembledSystem=FEMAssembly("EigenMode",2,1);
    if(frequency.NF==1)
        switch boundary.Axis
             case "x" 
                [MatA,MatB]=TFLF_Assembly_xBoundary(toolboxModel,boundary,Frequency);AssembledSystem.Matrix_A=MatA;AssembledSystem.Matrix_B=MatB;
            case "y"
                [MatA,MatB]=TFLF_Assembly_yBoundary(toolboxModel,boundary,Frequency);AssembledSystem.Matrix_A=MatA;AssembledSystem.Matrix_B=MatB;
            case "z"
                [MatA,MatB]=TFLF_Assembly_zBoundary(toolboxModel,boundary,Frequency);AssembledSystem.Matrix_A=MatA;AssembledSystem.Matrix_B=MatB;
        end
    else
        switch boundary.Axis
             case "x" 
                 for ii=1:frequency.NF,freq=frequency.Frequency(ii);[MatA,MatB]=TFLF_Assembly_xBoundary(toolboxModel,boundary,freq);AssembledSystem.Matrix_A{ii}=MatA;AssembledSystem.Matrix_B{ii}=MatB;end
            case "y"
                for ii=1:frequency.NF,freq=frequency.Frequency(ii);[MatA,MatB]=TFLF_Assembly_yBoundary(toolboxModel,boundary,freq);AssembledSystem.Matrix_A{ii}=MatA;AssembledSystem.Matrix_B{ii}=MatB;end
            case "z"
                for ii=1:frequency.NF,freq=frequency.Frequency(ii);[MatA,MatB]=TFLF_Assembly_zBoundary(toolboxModel,boundary,freq);AssembledSystem.Matrix_A{ii}=MatA;AssembledSystem.Matrix_B{ii}=MatB;end
        end
    end
end
%---------------------- x Propagation Assembly ----------------------------
function [MatA,MatB] = TFLF_Assembly_xBoundary(toolboxModel,boundary,Frequency),omega=2*pi*Frequency;
    II = zeros(120*9*numel(boundary.Facets),1); JJ = zeros(120*9*numel(boundary.Facets),1);SA = zeros(120*9*numel(boundary.Facets),1); SB =zeros(120*9*numel(boundary.Facets),1);counter=0;
    for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
        switch medium.Type
            case "Iso"
                [II,JJ,SA,SB,counter]=IsoAssembly_x(II,JJ,SA,SB,counter,toolboxModel,omega,element);
            case "Anis"
                [II,JJ,SA,SB,counter]=AnisAssembly_x(II,JJ,SA,SB,counter,toolboxModel,omega,element);
        end
    end
    NonZeros =nnz(II); SA=SA(1:NonZeros); II=II(1:NonZeros); JJ=JJ(1:NonZeros);SB=SB(1:NonZeros);
    MatA=sparse(II,JJ,SA);MatB=sparse(II,JJ,SB);
end
function [II,JJ,SA,SB,counter]=IsoAssembly_x(II,JJ,SA,SB,counter,toolboxModel,omega,element),GaussianQuadratture2D;medium=element.Medium2D;ElectromagneticConstants;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];z=[vertices.Z];De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;b(2)=(z(3)-z(1))/De;b(3)=(z(1)-z(2))/De;c(1)=(y(3)-y(2))/De;c(2)=(y(1)-y(3))/De;c(3)=(y(2)-y(1))/De;   
    Sh=zeros(3,3);Se=zeros(3,3);H=zeros(3,3);Th=zeros(3,3);Te=zeros(3,3);P=zeros(3,3);Tzh=zeros(3,3);Tze=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwy(1)=simp(1)*b(2)-simp(2)*b(1);wwy(2)=simp(2)*b(3)-simp(3)*b(2);wwy(3)=simp(3)*b(1)-simp(1)*b(3);
        wwz(1)=simp(1)*c(2)-simp(2)*c(1);wwz(2)=simp(2)*c(3)-simp(3)*c(2);wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wxy(1)=wwz(1);wxz(1)=-wwy(1);wxy(2)=wwz(2);wxz(2)=-wwy(2);wxy(3)=wwz(3);wxz(3)=-wwy(3);
        dxy(1)=c(1);dxz(1)=-b(1);dxy(2)=c(2);dxz(2)=-b(2);dxy(3)=c(3);dxz(3)=-b(3);
        wwy=wwy.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        for ii=1:3,for jj=1:3,Sh(ii,jj)=Sh(ii,jj)+Weights2D(kt)*(rw(ii)*(1/medium.Epsilon)*rw(jj));
                              Se(ii,jj)=Se(ii,jj)+Weights2D(kt)*(rw(ii)*(1/medium.Mu)*rw(jj));
                              H(ii,jj)=H(ii,jj)+Weights2D(kt)*(wwy(ii)*wxy(jj)+wwz(ii)*wxz(jj));
                              Th(ii,jj)=Th(ii,jj)+Weights2D(kt)*medium.Mu*(wwy(ii)*wwy(jj)+wwz(ii)*wwz(jj));
                              Te(ii,jj)=Te(ii,jj)+Weights2D(kt)*medium.Epsilon*(wwy(ii)*wwy(jj)+wwz(ii)*wwz(jj));
                              P(ii,jj)=P(ii,jj)+Weights2D(kt)*(simp(ii)*rw(jj));
                              Tze(ii,jj)=Tze(ii,jj)+Weights2D(kt)*(simp(ii)*medium.Epsilon*simp(jj));
                              Tzh(ii,jj)=Tzh(ii,jj)+Weights2D(kt)*(simp(ii)*medium.Mu*simp(jj));
                   end
        end
    end,Sh=Sh*Ae;Se=Se*Ae;H=H*Ae;Th=Th*Ae;Te=Te*Ae;P=P*Ae;Tzh=Tzh*Ae;Tze=Tzh*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;vi_E=vertices(ii).IndexE;vi_H=vertices(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;vj_E=vertices(jj).IndexE;vj_H=vertices(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=-omega*m0*si*sj*Th(ii,jj);
                                   counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=(1/(omega*e0))*si*sj*Sh(ii,jj);
            end
            if(ei_E~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=omega*e0*si*sj*Te(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-(1/(omega*m0))*si*sj*Se(ii,jj);
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=k0*si*sj*H(ii,jj);end
            if(ei_H~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_E;SB(counter)=k0*si*sj*H(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=vj_E;SA(counter)=-1i*omega*e0*Tze(ii,jj);end
            if(vi_E~=0 && ej_H~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=ej_H;SA(counter)=sj*P(ii,jj);end
            if(vi_H~=0 && vj_H~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=vj_H;SA(counter)=1i*omega*m0*Tzh(ii,jj);end
            if(vi_H~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=ej_E;SA(counter)=sj*P(ii,jj);end
        end
    end
end
function [II,JJ,SA,SB,counter]=AnisAssembly_x(II,JJ,SA,SB,counter,toolboxModel,omega,element),GaussianQuadratture2D;medium=element.Medium2D;ElectromagneticConstants;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    y=[vertices.Y];z=[vertices.Z];De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;b(2)=(z(3)-z(1))/De;b(3)=(z(1)-z(2))/De;c(1)=(y(3)-y(2))/De;c(2)=(y(1)-y(3))/De;c(3)=(y(2)-y(1))/De;   
    Sh=zeros(3,3);Se=zeros(3,3);H=zeros(3,3);Th=zeros(3,3);Te=zeros(3,3);P=zeros(3,3);Tzh=zeros(3,3);Tze=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwy(1)=simp(1)*b(2)-simp(2)*b(1);wwy(2)=simp(2)*b(3)-simp(3)*b(2);wwy(3)=simp(3)*b(1)-simp(1)*b(3);
        wwz(1)=simp(1)*c(2)-simp(2)*c(1);wwz(2)=simp(2)*c(3)-simp(3)*c(2);wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wxy(1)=wwz(1);wxz(1)=-wwy(1);wxy(2)=wwz(2);wxz(2)=-wwy(2);wxy(3)=wwz(3);wxz(3)=-wwy(3);
        dxy(1)=c(1);dxz(1)=-b(1);dxy(2)=c(2);dxz(2)=-b(2);dxy(3)=c(3);dxz(3)=-b(3);
        wwy=wwy.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        for ii=1:3,for jj=1:3,Sh(ii,jj)=Sh(ii,jj)+Weights2D(kt)*(rw(ii)*(1/medium.Epsilon(1,1))*rw(jj));
                              Se(ii,jj)=Se(ii,jj)+Weights2D(kt)*(rw(ii)*(1/medium.Mu(1,1))*rw(jj));
                              H(ii,jj)=H(ii,jj)+Weights2D(kt)*(wwy(ii)*wxy(jj)+wwz(ii)*wxz(jj));
                              Th(ii,jj)=Th(ii,jj)+Weights2D(kt)*(wwy(ii)*medium.Mu(2,2)*wwy(jj)+wwy(ii)*medium.Mu(2,3)*wwz(jj)+wwz(ii)*medium.Mu(3,3)*wwz(jj)+wwz(ii)*medium.Mu(3,2)*wwy(jj));
                              Te(ii,jj)=Te(ii,jj)+Weights2D(kt)*(wwy(ii)*medium.Epsilon(2,2)*wwy(jj)+wwy(ii)*medium.Epsilon(2,3)*wwz(jj)+wwz(ii)*medium.Epsilon(3,3)*wwz(jj)+wwz(ii)*medium.Epsilon(3,2)*wwy(jj));
                              P(ii,jj)=P(ii,jj)+Weights2D(kt)*(simp(ii)*rw(jj));
                              Tze(ii,jj)=Tze(ii,jj)+Weights2D(kt)*(simp(ii)*medium.Epsilon(1,1)*simp(jj));
                              Tzh(ii,jj)=Tzh(ii,jj)+Weights2D(kt)*(simp(ii)*medium.Mu(1,1)*simp(jj));
                   end
        end
    end,Sh=Sh*Ae;Se=Se*Ae;H=H*Ae;Th=Th*Ae;Te=Te*Ae;P=P*Ae;Tzh=Tzh*Ae;Tze=Tzh*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;vi_E=vertices(ii).IndexE;vi_H=vertices(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;vj_E=vertices(jj).IndexE;vj_H=vertices(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=-omega*m0*si*sj*Th(ii,jj);
                                   counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=(1/(omega*e0))*si*sj*Sh(ii,jj);
            end
            if(ei_E~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=omega*e0*si*sj*Te(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-(1/(omega*m0))*si*sj*Se(ii,jj);
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=k0*si*sj*H(ii,jj);end
            if(ei_H~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_E;SB(counter)=k0*si*sj*H(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=vj_E;SA(counter)=-1i*omega*e0*Tze(ii,jj);end
            if(vi_E~=0 && ej_H~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=ej_H;SA(counter)=sj*P(ii,jj);end
            if(vi_H~=0 && vj_H~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=vj_H;SA(counter)=1i*omega*m0*Tzh(ii,jj);end
            if(vi_H~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=ej_E;SA(counter)=sj*P(ii,jj);end
        end
    end
end
%---------------------- y Propagation Assembly ----------------------------
function [MatA,MatB] = TFLF_Assembly_yBoundary(toolboxModel,boundary,Frequency),omega=2*pi*Frequency;
    II = zeros(120*9*numel(boundary.Facets),1); JJ = zeros(120*9*numel(boundary.Facets),1);SA = zeros(120*9*numel(boundary.Facets),1); SB =zeros(120*9*numel(boundary.Facets),1);counter=0;
    for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
        switch medium.Type
            case "Iso"
                [II,JJ,SA,SB,counter]=IsoAssembly_y(II,JJ,SA,SB,counter,toolboxModel,omega,element);
            case "Anis"
                [II,JJ,SA,SB,counter]=AnisAssembly_y(II,JJ,SA,SB,counter,toolboxModel,omega,element);
        end
    end
    NonZeros =nnz(II); SA=SA(1:NonZeros); II=II(1:NonZeros); JJ=JJ(1:NonZeros);SB=SB(1:NonZeros);
    MatA=sparse(II,JJ,SA);MatB=sparse(II,JJ,SB);
end
function [II,JJ,SA,SB,counter]=IsoAssembly_y(II,JJ,SA,SB,counter,toolboxModel,omega,element),GaussianQuadratture2D;medium=element.Medium2D;ElectromagneticConstants;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];z=[vertices.Z];De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;b(2)=(z(3)-z(1))/De;b(3)=(z(1)-z(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;   
    Sh=zeros(3,3);Se=zeros(3,3);H=zeros(3,3);Th=zeros(3,3);Te=zeros(3,3);P=zeros(3,3);Tzh=zeros(3,3);Tze=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwz(1)=simp(1)*c(2)-simp(2)*c(1);wwz(2)=simp(2)*c(3)-simp(3)*c(2);wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wyx(1)=wwz(1);wyz(1)=-wwx(1);wyx(2)=wwz(2);wyz(2)=-wwx(2);wyx(3)=wwz(3);wyz(3)=-wwx(3);
        dyx(1)=c(1);dyz(1)=-b(1);dyz(2)=c(2);dyz(2)=-b(2);dyz(3)=c(3);dyz(3)=-b(3);
        wwy=wwy.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        for ii=1:3,for jj=1:3,Sh(ii,jj)=Sh(ii,jj)+Weights2D(kt)*(rw(ii)*(1/medium.Epsilon)*rw(jj));
                              Se(ii,jj)=Se(ii,jj)+Weights2D(kt)*(rw(ii)*(1/medium.Mu)*rw(jj));
                              H(ii,jj)=H(ii,jj)+Weights2D(kt)*(wwx(ii)*wyx(jj)+wwz(ii)*wyz(jj));
                              Th(ii,jj)=Th(ii,jj)+Weights2D(kt)*medium.Mu*(wwx(ii)*wwx(jj)+wwz(ii)*wwz(jj));
                              Te(ii,jj)=Te(ii,jj)+Weights2D(kt)*medium.Epsilon*(wwx(ii)*wwx(jj)+wwz(ii)*wwz(jj));
                              P(ii,jj)=P(ii,jj)+Weights2D(kt)*(simp(ii)*rw(jj));
                              Tze(ii,jj)=Tze(ii,jj)+Weights2D(kt)*(simp(ii)*medium.Epsilon*simp(jj));
                              Tzh(ii,jj)=Tzh(ii,jj)+Weights2D(kt)*(simp(ii)*medium.Mu*simp(jj));
                   end
        end
    end,Sh=Sh*Ae;Se=Se*Ae;H=H*Ae;Th=Th*Ae;Te=Te*Ae;P=P*Ae;Tzh=Tzh*Ae;Tze=Tzh*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;vi_E=vertices(ii).IndexE;vi_H=vertices(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;vj_E=vertices(jj).IndexE;vj_H=vertices(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=-omega*m0*si*sj*Th(ii,jj);
                                   counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=(1/(omega*e0))*si*sj*Sh(ii,jj);
            end
            if(ei_E~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=omega*e0*si*sj*Te(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-(1/(omega*m0))*si*sj*Se(ii,jj);
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=k0*si*sj*H(ii,jj);end
            if(ei_H~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_E;SB(counter)=k0*si*sj*H(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=vj_E;SA(counter)=-1i*omega*e0*Tze(ii,jj);end
            if(vi_E~=0 && ej_H~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=ej_H;SA(counter)=sj*P(ii,jj);end
            if(vi_H~=0 && vj_H~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=vj_H;SA(counter)=1i*omega*m0*Tzh(ii,jj);end
            if(vi_H~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=ej_E;SA(counter)=sj*P(ii,jj);end
        end
    end
end
function [II,JJ,SA,SB,counter]=AnisAssembly_y(II,JJ,SA,SB,counter,toolboxModel,omega,element),GaussianQuadratture2D;medium=element.Medium2D;ElectromagneticConstants;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];z=[vertices.Z];De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');Ae=abs(De)/2;
    b(1)=(z(2)-z(3))/De;b(2)=(z(3)-z(1))/De;b(3)=(z(1)-z(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;   
    Sh=zeros(3,3);Se=zeros(3,3);H=zeros(3,3);Th=zeros(3,3);Te=zeros(3,3);P=zeros(3,3);Tzh=zeros(3,3);Tze=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwz(1)=simp(1)*c(2)-simp(2)*c(1);wwz(2)=simp(2)*c(3)-simp(3)*c(2);wwz(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wyx(1)=wwz(1);wyz(1)=-wwx(1);wyx(2)=wwz(2);wyz(2)=-wwx(2);wyx(3)=wwz(3);wyz(3)=-wwx(3);
        dyx(1)=c(1);dyz(1)=-b(1);dyz(2)=c(2);dyz(2)=-b(2);dyz(3)=c(3);dyz(3)=-b(3);
        wwx=wwx.*edgeLength;wwz=wwz.*edgeLength;rw=rw.*edgeLength;
        for ii=1:3,for jj=1:3,Sh(ii,jj)=Sh(ii,jj)+Weights2D(kt)*(rw(ii)*(1/medium.Epsilon(2,2))*rw(jj));
                              Se(ii,jj)=Se(ii,jj)+Weights2D(kt)*(rw(ii)*(1/medium.Mu(2,2))*rw(jj));
                              H(ii,jj)=H(ii,jj)+Weights2D(kt)*(wwx(ii)*wyx(jj)+wwz(ii)*wyz(jj));
                              Th(ii,jj)=Th(ii,jj)+Weights2D(kt)*(wwx(ii)*medium.Mu(1,1)*wwx(jj)+wwx(ii)*medium.Mu(1,3)*wwz(jj)+wwz(ii)*medium.Mu(3,3)*wwz(jj)+wwz(ii)*medium.Mu(3,1)*wwx(jj));
                              Te(ii,jj)=Te(ii,jj)+Weights2D(kt)*(wwx(ii)*medium.Epsilon(1,1)*wwx(jj)+wwx(ii)*medium.Epsilon(1,3)*wwz(jj)+wwz(ii)*medium.Epsilon(3,3)*wwz(jj)+wwz(ii)*medium.Epsilon(3,1)*wwx(jj));
                              P(ii,jj)=P(ii,jj)+Weights2D(kt)*(simp(ii)*rw(jj));
                              Tze(ii,jj)=Tze(ii,jj)+Weights2D(kt)*(simp(ii)*medium.Epsilon(2,2)*simp(jj));
                              Tzh(ii,jj)=Tzh(ii,jj)+Weights2D(kt)*(simp(ii)*medium.Mu(2,2)*simp(jj));
                   end
        end
    end,Sh=Sh*Ae;Se=Se*Ae;H=H*Ae;Th=Th*Ae;Te=Te*Ae;P=P*Ae;Tzh=Tzh*Ae;Tze=Tzh*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;vi_E=vertices(ii).IndexE;vi_H=vertices(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;vj_E=vertices(jj).IndexE;vj_H=vertices(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=-omega*m0*si*sj*Th(ii,jj);
                                   counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=(1/(omega*e0))*si*sj*Sh(ii,jj);
            end
            if(ei_E~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=omega*e0*si*sj*Te(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-(1/(omega*m0))*si*sj*Se(ii,jj);
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=k0*si*sj*H(ii,jj);end
            if(ei_H~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_E;SB(counter)=k0*si*sj*H(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=vj_E;SA(counter)=-1i*omega*e0*Tze(ii,jj);end
            if(vi_E~=0 && ej_H~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=ej_H;SA(counter)=sj*P(ii,jj);end
            if(vi_H~=0 && vj_H~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=vj_H;SA(counter)=1i*omega*m0*Tzh(ii,jj);end
            if(vi_H~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=ej_E;SA(counter)=sj*P(ii,jj);end
        end
    end
end
%---------------------- z Propagation Assembly ----------------------------
function [MatA,MatB] = TFLF_Assembly_zBoundary(toolboxModel,boundary,Frequency),omega=2*pi*Frequency;
    II = zeros(120*9*numel(boundary.Facets),1); JJ = zeros(120*9*numel(boundary.Facets),1);SA = zeros(120*9*numel(boundary.Facets),1); SB =zeros(120*9*numel(boundary.Facets),1);counter=0;
    for ie=1:numel(boundary.Facets),element=toolboxModel.Facets(boundary.Facets(ie));medium=element.Medium2D;
        switch medium.Type
            case "Iso"
                [II,JJ,SA,SB,counter]=IsoAssembly_z(II,JJ,SA,SB,counter,toolboxModel,omega,element);
            case "Anis"
                [II,JJ,SA,SB,counter]=AnisAssembly_z(II,JJ,SA,SB,counter,toolboxModel,omega,element);
        end
    end
    NonZeros =nnz(II); SA=SA(1:NonZeros); II=II(1:NonZeros); JJ=JJ(1:NonZeros);SB=SB(1:NonZeros);
    MatA=sparse(II,JJ,SA);MatB=sparse(II,JJ,SB);
end
function [II,JJ,SA,SB,counter]=IsoAssembly_z(II,JJ,SA,SB,counter,toolboxModel,omega,element),GaussianQuadratture2D;medium=element.Medium2D;ElectromagneticConstants;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];y=[vertices.Y];De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');Ae=abs(De)/2;
    b(1)=(y(2)-y(3))/De;b(2)=(y(3)-y(1))/De;b(3)=(y(1)-y(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;   
    Sh=zeros(3,3);Se=zeros(3,3);H=zeros(3,3);Th=zeros(3,3);Te=zeros(3,3);P=zeros(3,3);Tzh=zeros(3,3);Tze=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwy(1)=simp(1)*c(2)-simp(2)*c(1);wwy(2)=simp(2)*c(3)-simp(3)*c(2);wwy(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wzx(1)=wwy(1);wzy(1)=-wwx(1);wzx(2)=wwy(2);wzy(2)=-wwx(2);wzx(3)=wwy(3);wzy(3)=-www(3);
        wwx=wwx.*edgeLength;wwy=wwy.*edgeLength;rw=rw.*edgeLength;
        for ii=1:3,for jj=1:3,Sh(ii,jj)=Sh(ii,jj)+Weights2D(kt)*(rw(ii)*(1/medium.Epsilon)*rw(jj));
                              Se(ii,jj)=Se(ii,jj)+Weights2D(kt)*(rw(ii)*(1/medium.Mu)*rw(jj));
                              H(ii,jj)=H(ii,jj)+Weights2D(kt)*(wwx(ii)*wzx(jj)+wwy(ii)*wzy(jj));
                              Th(ii,jj)=Th(ii,jj)+Weights2D(kt)*medium.Mu*(wwx(ii)*wwx(jj)+wwy(ii)*wwy(jj));
                              Te(ii,jj)=Te(ii,jj)+Weights2D(kt)*medium.Epsilon*(wwx(ii)*wwx(jj)+wwy(ii)*wwy(jj));
                              P(ii,jj)=P(ii,jj)+Weights2D(kt)*(simp(ii)*rw(jj));
                              Tze(ii,jj)=Tze(ii,jj)+Weights2D(kt)*(simp(ii)*medium.Epsilon*simp(jj));
                              Tzh(ii,jj)=Tzh(ii,jj)+Weights2D(kt)*(simp(ii)*medium.Mu*simp(jj));
                   end
        end
    end,Sh=Sh*Ae;Se=Se*Ae;H=H*Ae;Th=Th*Ae;Te=Te*Ae;P=P*Ae;Tzh=Tzh*Ae;Tze=Tzh*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;vi_E=vertices(ii).IndexE;vi_H=vertices(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;vj_E=vertices(jj).IndexE;vj_H=vertices(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=-omega*m0*si*sj*Th(ii,jj);
                                   counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=(1/(omega*e0))*si*sj*Sh(ii,jj);
            end
            if(ei_E~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=omega*e0*si*sj*Te(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-(1/(omega*m0))*si*sj*Se(ii,jj);
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=k0*si*sj*H(ii,jj);end
            if(ei_H~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_E;SB(counter)=k0*si*sj*H(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=vj_E;SA(counter)=-1i*omega*e0*Tze(ii,jj);end
            if(vi_E~=0 && ej_H~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=ej_H;SA(counter)=sj*P(ii,jj);end
            if(vi_H~=0 && vj_H~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=vj_H;SA(counter)=1i*omega*m0*Tzh(ii,jj);end
            if(vi_H~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=ej_E;SA(counter)=sj*P(ii,jj);end
        end
    end
end
function [II,JJ,SA,SB,counter]=AnisAssembly_z(II,JJ,SA,SB,counter,toolboxModel,omega,element),GaussianQuadratture2D;medium=element.Medium2D;ElectromagneticConstants;k0=omega/c0;
    vertices=[toolboxModel.Vertices(element.Vertices)];edges=[toolboxModel.Edges(element.Edges)];edgeLength=[edges.Length];
    x=[vertices.X];y=[vertices.Y];De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');Ae=abs(De)/2;
    b(1)=(y(2)-y(3))/De;b(2)=(y(3)-y(1))/De;b(3)=(y(1)-y(2))/De;c(1)=(x(3)-x(2))/De;c(2)=(x(1)-x(3))/De;c(3)=(x(2)-x(1))/De;   
    Sh=zeros(3,3);Se=zeros(3,3);H=zeros(3,3);Th=zeros(3,3);Te=zeros(3,3);P=zeros(3,3);Tzh=zeros(3,3);Tze=zeros(3,3);
    for kt=1:numel(Weights2D),simp(1)=r1_2D(kt);simp(2)=r2_2D(kt);simp(3)=r3_2D(kt);
        wwx(1)=simp(1)*b(2)-simp(2)*b(1);wwx(2)=simp(2)*b(3)-simp(3)*b(2);wwx(3)=simp(3)*b(1)-simp(1)*b(3);
        wwy(1)=simp(1)*c(2)-simp(2)*c(1);wwy(2)=simp(2)*c(3)-simp(3)*c(2);wwy(3)=simp(3)*c(1)-simp(1)*c(3);
        rw(1)=2*b(1)*c(2)-2*b(2)*c(1);rw(2)=2*b(2)*c(3)-2*b(3)*c(2);rw(3)=2*b(3)*c(1)-2*b(1)*c(3);
        wzx(1)=wwy(1);wzy(1)=-wwx(1);wzx(2)=wwy(2);wzy(2)=-wwx(2);wzx(3)=wwy(3);wzy(3)=-wwx(3);
        wwx=wwx.*edgeLength;wwy=wwy.*edgeLength;rw=rw.*edgeLength;
        for ii=1:3,for jj=1:3,Sh(ii,jj)=Sh(ii,jj)+Weights2D(kt)*(rw(ii)*(1/medium.Epsilon(3,3))*rw(jj));
                              Se(ii,jj)=Se(ii,jj)+Weights2D(kt)*(rw(ii)*(1/medium.Mu(3,3))*rw(jj));
                              H(ii,jj)=H(ii,jj)+Weights2D(kt)*(wwx(ii)*wzx(jj)+wwy(ii)*wzy(jj));
                              Th(ii,jj)=Th(ii,jj)+Weights2D(kt)*(wwx(ii)*medium.Mu(1,1)*wwx(jj)+wwx(ii)*medium.Mu(1,2)*wwy(jj)+wwy(ii)*medium.Mu(2,2)*wwy(jj)+wwy(ii)*medium.Mu(2,1)*wwx(jj));
                              Te(ii,jj)=Te(ii,jj)+Weights2D(kt)*(wwx(ii)*medium.Epsilon(1,1)*wwx(jj)+wwx(ii)*medium.Epsilon(1,2)*wwy(jj)+wwy(ii)*medium.Epsilon(2,2)*wwy(jj)+wwy(ii)*medium.Epsilon(2,1)*wwx(jj));
                              P(ii,jj)=P(ii,jj)+Weights2D(kt)*(simp(ii)*rw(jj));
                              Tze(ii,jj)=Tze(ii,jj)+Weights2D(kt)*(simp(ii)*medium.Epsilon(3,3)*simp(jj));
                              Tzh(ii,jj)=Tzh(ii,jj)+Weights2D(kt)*(simp(ii)*medium.Mu(3,3)*simp(jj));
                   end
        end
    end,Sh=Sh*Ae;Se=Se*Ae;H=H*Ae;Th=Th*Ae;Te=Te*Ae;P=P*Ae;Tzh=Tzh*Ae;Tze=Tzh*Ae;
    for ii=1:3,si=element.EdgeSigns(ii);ei_E=edges(ii).IndexE;ei_H=edges(ii).IndexH;vi_E=vertices(ii).IndexE;vi_H=vertices(ii).IndexH;
        for jj=1:3,sj=element.EdgeSigns(jj);ej_E=edges(jj).IndexE;ej_H=edges(jj).IndexH;vj_E=vertices(jj).IndexE;vj_H=vertices(jj).IndexH;
            if(ei_H~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=-omega*m0*si*sj*Th(ii,jj);
                                   counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_H;SA(counter)=(1/(omega*e0))*si*sj*Sh(ii,jj);
            end
            if(ei_E~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=omega*e0*si*sj*Te(ii,jj);
                                   counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_E;SA(counter)=-(1/(omega*m0))*si*sj*Se(ii,jj);
            end
            if(ei_E~=0 && ej_H~=0),counter=counter+1;II(counter)=ei_E;JJ(counter)=ej_H;SB(counter)=k0*si*sj*H(ii,jj);end
            if(ei_H~=0 && ej_E~=0),counter=counter+1;II(counter)=ei_H;JJ(counter)=ej_E;SB(counter)=k0*si*sj*H(ii,jj);end
            if(vi_E~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=vj_E;SA(counter)=-1i*omega*e0*Tze(ii,jj);end
            if(vi_E~=0 && ej_H~=0),counter=counter+1;II(counter)=vi_E;JJ(counter)=ej_H;SA(counter)=sj*P(ii,jj);end
            if(vi_H~=0 && vj_H~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=vj_H;SA(counter)=1i*omega*m0*Tzh(ii,jj);end
            if(vi_H~=0 && vj_E~=0),counter=counter+1;II(counter)=vi_H;JJ(counter)=ej_E;SA(counter)=sj*P(ii,jj);end
        end
    end
end
%}