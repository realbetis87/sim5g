function TModel =Assembly_ExcitationModule(TModel),AssembledSystem=TModel.Assembled;Frequency=TModel.Frequency;bian_flag=BianCheck(TModel);
%---------------------------- Single Frequency-----------------------------
    if(Frequency.NF==1)
        %--------------------- Dirichlet Boundary -------------------------
        if(TModel.Assembled.IsDir)
            %----------------- Vector Initializations ---------------------
            if(TModel.Assembled.Bian)
            else
                I=zeros(400*TModel.NumberOfElements,1);IA=I;JA=I;TE=I;TB=I;AA=I;FF=I;TS=I;TG=I;TBC=I;counterA=0;
                                                      IB=I;JB=I;B=I;counterB=0;
            end
            for id=1:numel(TModel.Domains),domain=TModel.Domains(id);medium=domain.Medium;
                switch medium.Type
                    case "Iso",[IA,JA,TE,TB,AA,FF,TS,TG,TBC,IB,JB,B,counterA,counterB] = Excitation_Dir_Iso(IA,JA,TE,TB,AA,FF,TS,TG,TBC,IB,JB,B,counterA,counterB,TModel,domain);
                    case "Anis",[IA,JA,TE,TB,AA,FF,TS,TG,TBC,IB,JB,B,counterA,counterB] = Excitation_Anis_Iso(IA,JA,TE,TB,AA,FF,TS,TG,TBC,IB,JB,B,counterA,counterB,TModel,domain);
                    case "Bian",[IA,JA,TE,TB,AA,FF,TS,TG,TBC,P,TA,TC,IB,JB,B,counterA,counterB] = Excitation_Bian_Iso(IA,JA,TE,TB,AA,FF,TS,TG,TBC,P,TA,TC,IB,JB,B,counterA,counterB,TModel,domain);
                end
            end
            NonZeros=nnz(IA);IA=IA(1:NonZeros);JA=JA(1:NonZeros);TE=TE(1:NonZeros);TB=TB(1:NonZeros);AA=AA(1:NonZeros);FF=FF(1:NonZeros);TS=TS(1:NonZeros);TG=TG(1:NonZeros);TBC=TBC(1:NonZeros);
            NonZerosB=nnz(IB);IB=IB(1:NonZerosB);JB=JB(1:NonZerosB);B=B(1:NonZerosB);freq=TModel.Frequency.Frequency;
            TE=sparse(IA,JA,TE);TB=sparse(IA,JA,TB);A=sparse(IA,JA,AA);F=sparse(IA,JA,FF);TS=sparse(IA,JA,TS);TG=sparse(IA,JA,TG);TBC=sparse(IA,JA,TBC);B=sparse(IB,JB,B);TModel.Assembled=TModel.Assembled.Excitation(TE,TB,A,F,TS,TG,TBC,B,freq);
        else
        end

%---------------------------- Multi Frequency-----------------------------
    else
    end
end
%=================================DIR Excitation ==========================
function [IA,JA,TE,TB,AA,FF,TS,TG,TBC,IB,JB,B,counterA,counterB] = Excitation_Dir_Iso(varargin),GaussianQuadratture3D;ElectromagneticConstants;
    if(nargin==16),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};IB=varargin{10};JB=varargin{11};B=varargin{12};
                   counterA=varargin{13};counterB=varargin{14};TModel=varargin{15};domain(16);
                   domain=varargin{18};medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;freq=TModel.Frequency.Frequency;
    elseif(nargin==17),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};IB=varargin{10};JB=varargin{11};B=varargin{12};
                   counterA=varargin{13};counterB=varargin{14};TModel=varargin{15};domain(16);FreqIndex=varargin{17};medium=domain.Medium;freq=TModel.Frequency.Frequency(FreqIndex);
                   if(medium.IsDispersive),epsilon=medium.Epsilon(FreqIndex);mu=medium.Mu(FreqIndex);imu=mu^-1;
                   else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;
                   end
    end,k0=2*pi*freq/c0;
    for ie=1:numel(domain.Elements),element=TModel.Elements(domain.Elements(ie));disp("element " + num2str(ie) +" out of " + num2str(numel(domain.Elements)));
        b=element.Bs;c=element.Cs;d=element.Ds;Ve=element.Volume;edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];
        edgeSigns=[element.EdgeSigns];facetSigns=[element.FacetSigns];edgeLengths=[edges.Length];facetSurfaces=[facets.Surface];
        Tbc=zeros(6,6);Ts=zeros(6,6);Tg=zeros(6,6);Te=zeros(6,6);Tb=zeros(4,4);Fm=zeros(4,6);Fw=zeros(4,6);Am=zeros(6,4);Aw=zeros(6,4);
        for kt=1:numel(Weights),zeta(1)=Points(1,kt);zeta(2)=Points(2,kt);zeta(3)=Points(3,kt);zeta(4)=Points(4,kt);
                 %---------------------- Edge Basis Functions -------------
                 wx(1)=zeta(1)*b(2)-zeta(2)*b(1); wx(2)=zeta(1)*b(3)-zeta(3)*b(1);wx(3)=zeta(1)*b(4)-zeta(4)*b(1);wx(4)=zeta(2)*b(3)-zeta(3)*b(2);wx(5)=zeta(2)*b(4)-zeta(4)*b(2);wx(6)=zeta(3)*b(4)-zeta(4)*b(3);
                 wy(1)=zeta(1)*c(2)-zeta(2)*c(1); wy(2)=zeta(1)*c(3)-zeta(3)*c(1);wy(3)=zeta(1)*c(4)-zeta(4)*c(1);wy(4)=zeta(2)*c(3)-zeta(3)*c(2);wy(5)=zeta(2)*c(4)-zeta(4)*c(2);wy(6)=zeta(3)*c(4)-zeta(4)*c(3);
                 wz(1)=zeta(1)*d(2)-zeta(2)*d(1); wz(2)=zeta(1)*d(3)-zeta(3)*d(1);wz(3)=zeta(1)*d(4)-zeta(4)*d(1);wz(4)=zeta(2)*d(3)-zeta(3)*d(2);wz(5)=zeta(2)*d(4)-zeta(4)*d(2);wz(6)=zeta(3)*d(4)-zeta(4)*d(3);
                 %------------------ Edge Basis Functions Rotations -------
                 rwx(1)=2*(c(1)*d(2)-d(1)*c(2));rwx(2)=2*(c(1)*d(3)-d(1)*c(3));rwx(3)=2*(c(1)*d(4)-d(1)*c(4));rwx(4)=2*(c(2)*d(3)-d(2)*c(3));rwx(5)=2*(c(2)*d(4)-d(2)*c(4));rwx(6)=2*(c(3)*d(4)-d(3)*c(4));
                 rwy(1)=2*(d(1)*b(2)-b(1)*d(2));rwy(2)=2*(d(1)*b(3)-b(1)*d(3));rwy(3)=2*(d(1)*b(4)-b(1)*d(4));rwy(4)=2*(d(2)*b(3)-b(2)*d(3));rwy(5)=2*(d(2)*b(4)-b(2)*d(4));rwy(6)=2*(d(3)*b(4)-b(3)*d(4));
                 rwz(1)=2*(b(1)*c(2)-c(1)*b(2));rwz(2)=2*(b(1)*c(3)-c(1)*b(3));rwz(3)=2*(b(1)*c(4)-c(1)*b(4));rwz(4)=2*(b(2)*c(3)-c(2)*b(3));rwz(5)=2*(b(2)*c(4)-c(2)*b(4));rwz(6)=2*(b(3)*c(4)-c(3)*b(4));
                 %---------------------- Facet Basis Functions ------------
                 wfx(1)=2*zeta(3)*(c(1)*d(2)-c(2)*d(1))+2*zeta(1)*(c(2)*d(3)-c(3)*d(2))+2*zeta(2)*(c(3)*d(1)-c(1)*d(3));wfy(1)=2*zeta(3)*(d(1)*b(2)-d(2)*b(1))+2*zeta(1)*(d(2)*b(3)-d(3)*b(2))+2*zeta(2)*(d(3)*b(1)-d(1)*b(3));wfz(1)=2*zeta(3)*(b(1)*c(2)-b(2)*c(1))+2*zeta(1)*(b(2)*c(3)-b(3)*c(2))+2*zeta(2)*(b(3)*c(1)-b(1)*c(3));
                 wfx(2)=2*zeta(3)*(c(2)*d(4)-c(4)*d(2))+2*zeta(2)*(c(4)*d(3)-c(3)*d(4))+2*zeta(4)*(c(3)*d(2)-c(2)*d(3));wfy(2)=2*zeta(3)*(d(2)*b(4)-d(4)*b(2))+2*zeta(2)*(d(4)*b(3)-d(3)*b(4))+2*zeta(4)*(d(3)*b(2)-d(2)*b(3));wfz(2)=2*zeta(3)*(b(2)*c(4)-b(4)*c(2))+2*zeta(2)*(b(4)*c(3)-b(3)*c(4))+2*zeta(4)*(b(3)*c(2)-b(2)*c(3));
                 wfx(3)=2*zeta(1)*(c(3)*d(4)-c(4)*d(3))+2*zeta(3)*(c(4)*d(1)-c(1)*d(4))+2*zeta(4)*(c(1)*d(3)-c(3)*d(1));wfy(3)=2*zeta(1)*(d(3)*b(4)-d(4)*b(3))+2*zeta(3)*(d(4)*b(1)-d(1)*b(4))+2*zeta(4)*(d(1)*b(3)-d(3)*b(1));wfz(3)=2*zeta(1)*(b(3)*c(4)-b(4)*c(3))+2*zeta(3)*(b(4)*c(1)-b(1)*c(4))+2*zeta(4)*(b(1)*c(3)-b(3)*c(1));
                 wfx(4)=2*zeta(1)*(c(4)*d(2)-c(2)*d(4))+2*zeta(4)*(c(2)*d(1)-c(1)*d(2))+2*zeta(2)*(c(1)*d(4)-c(4)*d(1));wfy(4)=2*zeta(1)*(d(4)*b(2)-d(2)*b(4))+2*zeta(4)*(d(2)*b(1)-d(1)*b(2))+2*zeta(2)*(d(1)*b(4)-d(4)*b(1));wfz(4)=2*zeta(1)*(b(4)*c(2)-b(2)*c(4))+2*zeta(4)*(b(2)*c(1)-b(1)*c(2))+2*zeta(2)*(b(1)*c(4)-b(4)*c(1));
                 %===================== Basis Function Scaling ============
                 switch AssembledSystem.E_Scaling
                     case 1,wx=wx.*edgeLengths;wy=wy.*edgeLengths;wz=wz.*edgeLengths;rwx=rwx.*edgeLengths;rwy=rwy.*edgeLengths;rwy=rwy.*edgeLengths;rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch AssembledSystem.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);end
                 end
                 %==================== Local Matrices =====================
                 for ii=1:6,for jj=1:6,Te(ii,jj)=Te(ii,jj)+Weights(kt)*epsilon*(wx(ii)*wx(jj)+wy(ii)*wy(jj)+wz(ii)*wz(jj));end,end
                 for ii=1:4,for jj=1:4,Tb(ii,jj)=Tb(ii,jj)+Weights(kt)*(wfx(ii)*wfx(jj)+wfy(ii)*wfy(jj)+wfz(ii)*wfz(jj));end,end
                 for ii=1:4,for jj=1:6,Fm(ii,jj)=Fm(ii,jj)+Weights(kt)*(wfx(ii)*rwx(jj)+wfy(ii)*rwy(jj)+wfz(ii)*rwz(jj));end,end
                 for ii=1:6,for jj=1:4,Am(ii,jj)=Am(ii,jj)+Weights(kt)*imu*(rwx(ii)*wfx(jj)+rwy(ii)*wfy(jj)+rwz(ii)*wfz(jj));end,end
        end,Te=Te*Ve;Tb=Tb*Ve;Am=Am*Ve;Fm=Fm*Ve;
        %------------------------------------------------------------------
        for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                if(boundary.Type=="ABC"),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=1i*sqrt(epsilon/mu)*Ts;
                   switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end 
                elseif(boundary.Type=="GRA"),Tg = CalculateIntegral(TModel,element,facets(ii),b,c,d,ii);
                    if(boundary.Dispersive),cond=boundary.Param(FreqIndex);else,cond=boundary.Param;end,Tg=Tg.*cond;
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                elseif(boundary.Type=="IBC"),Tbc = CalculateIntegral(TModel,element,facets(ii),b,c,d,ii);
                    if(boundary.Dispersive),cond=boundary.Param(FreqIndex);else,cond=boundary.Param;end,Tbc=Tbc.*cond;
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                end
            end
        end
        %------------------------------------------------------------------
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);ki=edges(ii).KnownIndex;
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*k0*Te(ii,jj);end
            end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);ki=edges(ii).KnownIndex;
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*Am(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);ki=facets(ii).KnownIndex;
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*Fm(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);ki=facets(ii).KnownIndex;
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*k0*Tb(ii,jj);end
            end
        end
    end
end
function [IA,JA,TE,TB,AA,FF,TS,TG,TBC,IB,JB,B,counterA,counterB] = Excitation_Anis_Iso(varargin),GaussianQuadratture3D;ElectromagneticConstants;
    if(nargin==16),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};IB=varargin{10};JB=varargin{11};B=varargin{12};
                   counterA=varargin{13};counterB=varargin{14};TModel=varargin{15};domain(16);
                   domain=varargin{18};medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;freq=TModel.Frequency.Frequency;
    elseif(nargin==17),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};IB=varargin{10};JB=varargin{11};B=varargin{12};
                   counterA=varargin{13};counterB=varargin{14};TModel=varargin{15};domain(16);FreqIndex=varargin{17};medium=domain.Medium;freq=TModel.Frequency.Frequency(FreqIndex);
                   if(medium.IsDispersive),epsilon=medium.Epsilon{FreqIndex};mu=medium.Mu{FreqIndex};imu=mu^-1;
                   else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;
                   end
    end,k0=2*pi*freq/c0;
    for ie=1:numel(domain.Elements),element=TModel.Elements(domain.Elements(ie));disp("element " + num2str(ie) +" out of " + num2str(numel(domain.Elements)));
        b=element.Bs;c=element.Cs;d=element.Ds;Ve=element.Volume;edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];
        edgeSigns=[element.EdgeSigns];facetSigns=[element.FacetSigns];edgeLengths=[edges.Length];facetSurfaces=[facets.Surface];
        Tbc=zeros(6,6);Ts=zeros(6,6);Tg=zeros(6,6);Te=zeros(6,6);Tb=zeros(4,4);Fm=zeros(4,6);Fw=zeros(4,6);Am=zeros(6,4);Aw=zeros(6,4);
        for kt=1:numel(Weights),zeta(1)=Points(1,kt);zeta(2)=Points(2,kt);zeta(3)=Points(3,kt);zeta(4)=Points(4,kt);
                 %---------------------- Edge Basis Functions -------------
                 wx(1)=zeta(1)*b(2)-zeta(2)*b(1); wx(2)=zeta(1)*b(3)-zeta(3)*b(1);wx(3)=zeta(1)*b(4)-zeta(4)*b(1);wx(4)=zeta(2)*b(3)-zeta(3)*b(2);wx(5)=zeta(2)*b(4)-zeta(4)*b(2);wx(6)=zeta(3)*b(4)-zeta(4)*b(3);
                 wy(1)=zeta(1)*c(2)-zeta(2)*c(1); wy(2)=zeta(1)*c(3)-zeta(3)*c(1);wy(3)=zeta(1)*c(4)-zeta(4)*c(1);wy(4)=zeta(2)*c(3)-zeta(3)*c(2);wy(5)=zeta(2)*c(4)-zeta(4)*c(2);wy(6)=zeta(3)*c(4)-zeta(4)*c(3);
                 wz(1)=zeta(1)*d(2)-zeta(2)*d(1); wz(2)=zeta(1)*d(3)-zeta(3)*d(1);wz(3)=zeta(1)*d(4)-zeta(4)*d(1);wz(4)=zeta(2)*d(3)-zeta(3)*d(2);wz(5)=zeta(2)*d(4)-zeta(4)*d(2);wz(6)=zeta(3)*d(4)-zeta(4)*d(3);
                 %------------------ Edge Basis Functions Rotations -------
                 rwx(1)=2*(c(1)*d(2)-d(1)*c(2));rwx(2)=2*(c(1)*d(3)-d(1)*c(3));rwx(3)=2*(c(1)*d(4)-d(1)*c(4));rwx(4)=2*(c(2)*d(3)-d(2)*c(3));rwx(5)=2*(c(2)*d(4)-d(2)*c(4));rwx(6)=2*(c(3)*d(4)-d(3)*c(4));
                 rwy(1)=2*(d(1)*b(2)-b(1)*d(2));rwy(2)=2*(d(1)*b(3)-b(1)*d(3));rwy(3)=2*(d(1)*b(4)-b(1)*d(4));rwy(4)=2*(d(2)*b(3)-b(2)*d(3));rwy(5)=2*(d(2)*b(4)-b(2)*d(4));rwy(6)=2*(d(3)*b(4)-b(3)*d(4));
                 rwz(1)=2*(b(1)*c(2)-c(1)*b(2));rwz(2)=2*(b(1)*c(3)-c(1)*b(3));rwz(3)=2*(b(1)*c(4)-c(1)*b(4));rwz(4)=2*(b(2)*c(3)-c(2)*b(3));rwz(5)=2*(b(2)*c(4)-c(2)*b(4));rwz(6)=2*(b(3)*c(4)-c(3)*b(4));
                 %---------------------- Facet Basis Functions ------------
                 wfx(1)=2*zeta(3)*(c(1)*d(2)-c(2)*d(1))+2*zeta(1)*(c(2)*d(3)-c(3)*d(2))+2*zeta(2)*(c(3)*d(1)-c(1)*d(3));wfy(1)=2*zeta(3)*(d(1)*b(2)-d(2)*b(1))+2*zeta(1)*(d(2)*b(3)-d(3)*b(2))+2*zeta(2)*(d(3)*b(1)-d(1)*b(3));wfz(1)=2*zeta(3)*(b(1)*c(2)-b(2)*c(1))+2*zeta(1)*(b(2)*c(3)-b(3)*c(2))+2*zeta(2)*(b(3)*c(1)-b(1)*c(3));
                 wfx(2)=2*zeta(3)*(c(2)*d(4)-c(4)*d(2))+2*zeta(2)*(c(4)*d(3)-c(3)*d(4))+2*zeta(4)*(c(3)*d(2)-c(2)*d(3));wfy(2)=2*zeta(3)*(d(2)*b(4)-d(4)*b(2))+2*zeta(2)*(d(4)*b(3)-d(3)*b(4))+2*zeta(4)*(d(3)*b(2)-d(2)*b(3));wfz(2)=2*zeta(3)*(b(2)*c(4)-b(4)*c(2))+2*zeta(2)*(b(4)*c(3)-b(3)*c(4))+2*zeta(4)*(b(3)*c(2)-b(2)*c(3));
                 wfx(3)=2*zeta(1)*(c(3)*d(4)-c(4)*d(3))+2*zeta(3)*(c(4)*d(1)-c(1)*d(4))+2*zeta(4)*(c(1)*d(3)-c(3)*d(1));wfy(3)=2*zeta(1)*(d(3)*b(4)-d(4)*b(3))+2*zeta(3)*(d(4)*b(1)-d(1)*b(4))+2*zeta(4)*(d(1)*b(3)-d(3)*b(1));wfz(3)=2*zeta(1)*(b(3)*c(4)-b(4)*c(3))+2*zeta(3)*(b(4)*c(1)-b(1)*c(4))+2*zeta(4)*(b(1)*c(3)-b(3)*c(1));
                 wfx(4)=2*zeta(1)*(c(4)*d(2)-c(2)*d(4))+2*zeta(4)*(c(2)*d(1)-c(1)*d(2))+2*zeta(2)*(c(1)*d(4)-c(4)*d(1));wfy(4)=2*zeta(1)*(d(4)*b(2)-d(2)*b(4))+2*zeta(4)*(d(2)*b(1)-d(1)*b(2))+2*zeta(2)*(d(1)*b(4)-d(4)*b(1));wfz(4)=2*zeta(1)*(b(4)*c(2)-b(2)*c(4))+2*zeta(4)*(b(2)*c(1)-b(1)*c(2))+2*zeta(2)*(b(1)*c(4)-b(4)*c(1));
                 %===================== Basis Function Scaling ============
                 switch AssembledSystem.E_Scaling
                     case 1,wx=wx.*edgeLengths;wy=wy.*edgeLengths;wz=wz.*edgeLengths;rwx=rwx.*edgeLengths;rwy=rwy.*edgeLengths;rwy=rwy.*edgeLengths;rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch AssembledSystem.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);end
                 end
                 %------------------------ (mr^-1)*wf----------------------
                 mwfx(1)=imu(1,1)*wfx(1)+imu(1,2)*wfy(1)+imu(1,3)*wfz(1);mwfy(1)=imu(2,1)*wfx(1)+imu(2,2)*wfy(1)+imu(2,3)*wfz(1);mwfz(1)=imu(3,1)*wfx(1)+imu(3,2)*wfy(1)+imu(3,3)*wfz(1);
                 mwfx(2)=imu(1,1)*wfx(2)+imu(1,2)*wfy(2)+imu(1,3)*wfz(2);mwfy(2)=imu(2,1)*wfx(2)+imu(2,2)*wfy(2)+imu(2,3)*wfz(2);mwfz(2)=imu(3,1)*wfx(2)+imu(3,2)*wfy(2)+imu(3,3)*wfz(2);
                 mwfx(3)=imu(1,1)*wfx(3)+imu(1,2)*wfy(3)+imu(1,3)*wfz(3);mwfy(3)=imu(2,1)*wfx(3)+imu(2,2)*wfy(3)+imu(2,3)*wfz(3);mwfz(3)=imu(3,1)*wfx(3)+imu(3,2)*wfy(3)+imu(3,3)*wfz(3);
                 mwfx(4)=imu(1,1)*wfx(4)+imu(1,2)*wfy(4)+imu(1,3)*wfz(4);mwfy(4)=imu(2,1)*wfx(4)+imu(2,2)*wfy(4)+imu(2,3)*wfz(4);mwfz(4)=imu(3,1)*wfx(4)+imu(3,2)*wfy(4)+imu(3,3)*wfz(4);
                 %-------------- x cross Edge Basis Functions -------------
                 xwx(1)=0;xwy(1)=-wz(1);xwz(1)=wy(1);xwx(2)=0;xwy(2)=-wz(2);xwz(2)=wy(2);xwx(3)=0;xwy(3)=-wz(3);xwz(3)=wy(3);
                 xwx(4)=0;xwy(4)=-wz(4);xwz(4)=wy(4);xwx(5)=0;xwy(5)=-wz(5);xwz(5)=wy(5);xwx(6)=0;xwy(6)=-wz(6);xwz(6)=wy(6);
                 %-------------- x cross mr^-1 Facet Basis Functions ------
                 xwfx(1)=0;xwfy(1)=-mwfz(1);xwfz(1)=mwfy(1);xwfx(2)=0;xwfy(2)=-mwfz(2);xwfz(2)=mwfy(2);
                 xwfx(3)=0;xwfy(3)=-mwfz(3);xwfz(3)=mwfy(3);xwfx(4)=0;xwfy(4)=-mwfz(4);xwfz(4)=mwfy(4);
                 %==================== Local Matrices =====================
                 for ii=1:6
                     for jj=1:6,Te(ii,jj)=Te(ii,jj)+Weights(kt)*(wx(ii)*epsilon(1,1)*wx(jj) +wx(ii)*epsilon(1,2)*wy(jj) +wx(ii)*epsilon(1,3)*wz(jj)...
                                                                  +wy(ii)*epsilon(2,1)*wx(jj)+wy(ii)*epsilon(2,2)*wy(jj)+wy(ii)*epsilon(2,3)*wz(jj)...
                                                                  +wz(ii)*epsilon(3,1)*wx(jj)+wz(ii)*epsilon(3,2)*wy(jj)+wz(ii)*epsilon(3,3)*wz(jj));
                     end
                 end
                 for ii=1:4,for jj=1:4,Tb(ii,jj)=Tb(ii,jj)+Weights(kt)*(wfx(ii)*wfx(jj)+wfy(ii)*wfy(jj)+wfz(ii)*wfz(jj));end,end
                 for ii=1:4,for jj=1:6,Fm(ii,jj)=Fm(ii,jj)+Weights(kt)*(wfx(ii)*rwx(jj)+wfy(ii)*rwy(jj)+wfz(ii)*rwz(jj));end,end
                 for ii=1:6,for jj=1:4,Am(ii,jj)=Am(ii,jj)+Weights(kt)*(rwx(ii)*mwfx(jj)+rwy(ii)*mwfy(jj)+rwz(ii)*mwfz(jj));end,end
        end,Te=Te*Ve;Tb=Tb*Ve;Am=Am*Ve;Fm=Fm*Ve;Fw=Fw*Ve;Aw=Aw*Ve;
        %------------------------------------------------------------------
        for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                if(boundary.Type=="ABC"),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=1i*sqrt(epsilon/mu)*Ts;
                   switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end 
                elseif(boundary.Type=="GRA"),Tg = CalculateIntegral(TModel,element,facets(ii),b,c,d,ii);
                    if(boundary.Dispersive),cond=boundary.Param(FreqIndex);else,cond=boundary.Param;end,Tg=Tg.*cond;
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                elseif(boundary.Type=="IBC"),Tbc = CalculateIntegral(TModel,element,facets(ii),b,c,d,ii);
                    if(boundary.Dispersive),cond=boundary.Param(FreqIndex);else,cond=boundary.Param;end,Tbc=Tbc.*cond;
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                end
            end
        end
        %------------------------------------------------------------------
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);ki=edges(ii).KnownIndex;
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*k0*Te(ii,jj);end
            end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);ki=edges(ii).KnownIndex;
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*Am(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);ki=facets(ii).KnownIndex;
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*Fm(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);ki=facets(ii).KnownIndex;
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*k0*Tb(ii,jj);end
            end
        end
    end
end
function [IA,JA,TE,TB,AA,FF,TS,TG,TBC,P,TA,TC,IB,JB,B,counterA,counterB] = Excitation_Bian_Iso(varargin),GaussianQuadratture3D;ElectromagneticConstants;
    if(nargin==19),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
               TS=varargin{7};TG=varargin{8};TBC=varargin{9};P=varargin{10};TA=varargin{11};TC=varargin{12};IB=varargin{13};
               JB=varargin{14};B=varargin{15};counterA=varargin{16};counterB=varargin{17};TModel.varargin{18};domain=varargin{19};
               medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;ksi=medium.Ksi;zita=medium.Zita;freq=TModel.Frequency.Frequency;
    elseif(nargin==20),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
               TS=varargin{7};TG=varargin{8};TBC=varargin{9};P=varargin{10};TA=varargin{11};TC=varargin{12};IB=varargin{13};
               JB=varargin{14};B=varargin{15};counterA=varargin{16};counterB=varargin{17};TModel.varargin{18};domain=varargin{19};FreqIndex=varargin{20};freq=TModel.Frequency.Frequency(FreqIndex);
               medium=domain.Medium;
               if(medium.IsDispersive),epsilon=medium.Epsilon{FreqIndex};mu=medium.Mu{FreqIndex};imu=mu^-1;ksi=medium.Ksi{FreqIndex};zita=medium.Zita{FreqIndex};
               else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;ksi=medium.Ksi;zita=medium.Zita;
               end
    end,k0=2*pi*freq/c0;
    for ie=1:numel(domain.Elements),element=TModel.Elements(domain.Elements(ie));
        b=element.Bs;c=element.Cs;d=element.Ds;Ve=element.Volume;edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];
        edgeSigns=[element.EdgeSigns];facetSigns=[element.FacetSigns];edgeLengths=[edges.Length];facetSurfaces=[facets.Surface];
        Tbc=zeros(6,6);Ts=zeros(6,6);Tg=zeros(6,6);Te=zeros(6,6);Tb=zeros(4,4);Fm=zeros(4,6);Fw=zeros(4,6);Am=zeros(6,4);Aw=zeros(6,4);Pm=zeros(6,6);Km=zeros(6,6);Tam=zeros(6,6);Tcm=zeros(6,4);
        for kt=1:numel(Weights),zeta(1)=Points(1,kt);zeta(2)=Points(2,kt);zeta(3)=Points(3,kt);zeta(4)=Points(4,kt);
             %---------------------- Edge Basis Functions -------------
                 wx(1)=zeta(1)*b(2)-zeta(2)*b(1); wx(2)=zeta(1)*b(3)-zeta(3)*b(1);wx(3)=zeta(1)*b(4)-zeta(4)*b(1);wx(4)=zeta(2)*b(3)-zeta(3)*b(2);wx(5)=zeta(2)*b(4)-zeta(4)*b(2);wx(6)=zeta(3)*b(4)-zeta(4)*b(3);
                 wy(1)=zeta(1)*c(2)-zeta(2)*c(1); wy(2)=zeta(1)*c(3)-zeta(3)*c(1);wy(3)=zeta(1)*c(4)-zeta(4)*c(1);wy(4)=zeta(2)*c(3)-zeta(3)*c(2);wy(5)=zeta(2)*c(4)-zeta(4)*c(2);wy(6)=zeta(3)*c(4)-zeta(4)*c(3);
                 wz(1)=zeta(1)*d(2)-zeta(2)*d(1); wz(2)=zeta(1)*d(3)-zeta(3)*d(1);wz(3)=zeta(1)*d(4)-zeta(4)*d(1);wz(4)=zeta(2)*d(3)-zeta(3)*d(2);wz(5)=zeta(2)*d(4)-zeta(4)*d(2);wz(6)=zeta(3)*d(4)-zeta(4)*d(3);
                 %------------------ Edge Basis Functions Rotations -------
                 rwx(1)=2*(c(1)*d(2)-d(1)*c(2));rwx(2)=2*(c(1)*d(3)-d(1)*c(3));rwx(3)=2*(c(1)*d(4)-d(1)*c(4));rwx(4)=2*(c(2)*d(3)-d(2)*c(3));rwx(5)=2*(c(2)*d(4)-d(2)*c(4));rwx(6)=2*(c(3)*d(4)-d(3)*c(4));
                 rwy(1)=2*(d(1)*b(2)-b(1)*d(2));rwy(2)=2*(d(1)*b(3)-b(1)*d(3));rwy(3)=2*(d(1)*b(4)-b(1)*d(4));rwy(4)=2*(d(2)*b(3)-b(2)*d(3));rwy(5)=2*(d(2)*b(4)-b(2)*d(4));rwy(6)=2*(d(3)*b(4)-b(3)*d(4));
                 rwz(1)=2*(b(1)*c(2)-c(1)*b(2));rwz(2)=2*(b(1)*c(3)-c(1)*b(3));rwz(3)=2*(b(1)*c(4)-c(1)*b(4));rwz(4)=2*(b(2)*c(3)-c(2)*b(3));rwz(5)=2*(b(2)*c(4)-c(2)*b(4));rwz(6)=2*(b(3)*c(4)-c(3)*b(4));
                 %---------------------- Facet Basis Functions ------------
                 wfx(1)=2*zeta(3)*(c(1)*d(2)-c(2)*d(1))+2*zeta(1)*(c(2)*d(3)-c(3)*d(2))+2*zeta(2)*(c(3)*d(1)-c(1)*d(3));wfy(1)=2*zeta(3)*(d(1)*b(2)-d(2)*b(1))+2*zeta(1)*(d(2)*b(3)-d(3)*b(2))+2*zeta(2)*(d(3)*b(1)-d(1)*b(3));wfz(1)=2*zeta(3)*(b(1)*c(2)-b(2)*c(1))+2*zeta(1)*(b(2)*c(3)-b(3)*c(2))+2*zeta(2)*(b(3)*c(1)-b(1)*c(3));
                 wfx(2)=2*zeta(3)*(c(2)*d(4)-c(4)*d(2))+2*zeta(2)*(c(4)*d(3)-c(3)*d(4))+2*zeta(4)*(c(3)*d(2)-c(2)*d(3));wfy(2)=2*zeta(3)*(d(2)*b(4)-d(4)*b(2))+2*zeta(2)*(d(4)*b(3)-d(3)*b(4))+2*zeta(4)*(d(3)*b(2)-d(2)*b(3));wfz(2)=2*zeta(3)*(b(2)*c(4)-b(4)*c(2))+2*zeta(2)*(b(4)*c(3)-b(3)*c(4))+2*zeta(4)*(b(3)*c(2)-b(2)*c(3));
                 wfx(3)=2*zeta(1)*(c(3)*d(4)-c(4)*d(3))+2*zeta(3)*(c(4)*d(1)-c(1)*d(4))+2*zeta(4)*(c(1)*d(3)-c(3)*d(1));wfy(3)=2*zeta(1)*(d(3)*b(4)-d(4)*b(3))+2*zeta(3)*(d(4)*b(1)-d(1)*b(4))+2*zeta(4)*(d(1)*b(3)-d(3)*b(1));wfz(3)=2*zeta(1)*(b(3)*c(4)-b(4)*c(3))+2*zeta(3)*(b(4)*c(1)-b(1)*c(4))+2*zeta(4)*(b(1)*c(3)-b(3)*c(1));
                 wfx(4)=2*zeta(1)*(c(4)*d(2)-c(2)*d(4))+2*zeta(4)*(c(2)*d(1)-c(1)*d(2))+2*zeta(2)*(c(1)*d(4)-c(4)*d(1));wfy(4)=2*zeta(1)*(d(4)*b(2)-d(2)*b(4))+2*zeta(4)*(d(2)*b(1)-d(1)*b(2))+2*zeta(2)*(d(1)*b(4)-d(4)*b(1));wfz(4)=2*zeta(1)*(b(4)*c(2)-b(2)*c(4))+2*zeta(4)*(b(2)*c(1)-b(1)*c(2))+2*zeta(2)*(b(1)*c(4)-b(4)*c(1));
                 %===================== Basis Function Scaling ============
                 switch AssembledSystem.E_Scaling
                     case 1,wx=wx.*edgeLengths;wy=wy.*edgeLengths;wz=wz.*edgeLengths;rwx=rwx.*edgeLengths;rwy=rwy.*edgeLengths;rwy=rwy.*edgeLengths;rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch AssembledSystem.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);end
                 end
                 %-------------  (mr^-1) * Facet Basis Functions-----------
                 mwfx(1)=imu(1,1)*wfx(1)+imu(1,2)*wfy(1)+imu(1,3)*wfz(1);mwfy(1)=imu(2,1)*wfx(1)+imu(2,2)*wfy(1)+imu(2,3)*wfz(1);mwfz(1)=imu(3,1)*wfx(1)+imu(3,2)*wfy(1)+imu(3,3)*wfz(1);
                 mwfx(2)=imu(1,1)*wfx(2)+imu(1,2)*wfy(2)+imu(1,3)*wfz(2);mwfy(2)=imu(2,1)*wfx(2)+imu(2,2)*wfy(2)+imu(2,3)*wfz(2);mwfz(1)=imu(3,1)*wfx(2)+imu(3,2)*wfy(2)+imu(3,3)*wfz(2);
                 mwfx(3)=imu(1,1)*wfx(3)+imu(1,2)*wfy(3)+imu(1,3)*wfz(3);mwfy(3)=imu(2,1)*wfx(3)+imu(2,2)*wfy(3)+imu(2,3)*wfz(3);mwfz(1)=imu(3,1)*wfx(3)+imu(3,2)*wfy(3)+imu(3,3)*wfz(3);
                 mwfx(4)=imu(1,1)*wfx(4)+imu(1,2)*wfy(4)+imu(1,3)*wfz(4);mwfy(4)=imu(2,1)*wfx(4)+imu(2,2)*wfy(4)+imu(2,3)*wfz(4);mwfz(1)=imu(3,1)*wfx(4)+imu(3,2)*wfy(4)+imu(3,3)*wfz(4);
                %----------- (mr^-1)*zita* Edge Basis Functions -----------
                mzwx(1)=imz(1,1)*wx(1)+imz(1,2)*wy(1)+imz(1,3)*wz(1);mzwy(1)=imz(2,1)*wx(1)+imz(2,2)*wy(1)+imz(2,3)*wz(1);mzwz(1)=imz(3,1)*wx(1)+imz(3,2)*wy(1)+imz(3,3)*wz(1);
                mzwx(2)=imz(1,1)*wx(2)+imz(1,2)*wy(2)+imz(1,3)*wz(2);mzwy(2)=imz(2,1)*wx(2)+imz(2,2)*wy(2)+imz(2,3)*wz(2);mzwz(2)=imz(3,1)*wx(2)+imz(3,2)*wy(2)+imz(3,3)*wz(2);
                mzwx(3)=imz(1,1)*wx(3)+imz(1,2)*wy(3)+imz(1,3)*wz(3);mzwy(3)=imz(2,1)*wx(3)+imz(2,2)*wy(3)+imz(2,3)*wz(3);mzwz(3)=imz(3,1)*wx(3)+imz(3,2)*wy(3)+imz(3,3)*wz(3);
                mzwx(4)=imz(1,1)*wx(4)+imz(1,2)*wy(4)+imz(1,3)*wz(4);mzwy(4)=imz(2,1)*wx(4)+imz(2,2)*wy(4)+imz(2,3)*wz(4);mzwz(4)=imz(3,1)*wx(4)+imz(3,2)*wy(4)+imz(3,3)*wz(4);
                mzwx(5)=imz(1,1)*wx(5)+imz(1,2)*wy(5)+imz(1,3)*wz(5);mzwy(5)=imz(2,1)*wx(5)+imz(2,2)*wy(5)+imz(2,3)*wz(5);mzwz(5)=imz(3,1)*wx(5)+imz(3,2)*wy(5)+imz(3,3)*wz(5);
                mzwx(6)=imz(1,1)*wx(6)+imz(1,2)*wy(6)+imz(1,3)*wz(6);mzwy(6)=imz(2,1)*wx(6)+imz(2,2)*wy(6)+imz(2,3)*wz(6);mzwz(6)=imz(3,1)*wx(6)+imz(3,2)*wy(6)+imz(3,3)*wz(6);
                %------------- ksi *(mr^-1)* FacetBasis Functions ---------
                kmvx(1)=kim(1,1)*wfx(1)+kim(1,2)*wfy(1)+kim(1,3)*wfz(1);kmvy(1)=kim(2,1)*wfx(1)+kim(2,2)*wfy(1)+kim(2,3)*wfz(1);kmvz(1)=kim(3,1)*wfx(1)+kim(3,2)*wfy(1)+kim(3,3)*wfz(1);
                kmvx(2)=kim(1,1)*wfx(2)+kim(1,2)*wfy(2)+kim(1,3)*wfz(2);kmvy(2)=kim(2,1)*wfx(2)+kim(2,2)*wfy(2)+kim(2,3)*wfz(2);kmvz(2)=kim(3,1)*wfx(2)+kim(3,2)*wfy(2)+kim(3,3)*wfz(2);
                kmvx(3)=kim(1,1)*wfx(3)+kim(1,2)*wfy(3)+kim(1,3)*wfz(3);kmvy(3)=kim(2,1)*wfx(3)+kim(2,2)*wfy(3)+kim(2,3)*wfz(3);kmvz(3)=kim(3,1)*wfx(3)+kim(3,2)*wfy(3)+kim(3,3)*wfz(3);
                kmvx(4)=kim(1,1)*wfx(4)+kim(1,2)*wfy(4)+kim(1,3)*wfz(4);kmvy(4)=kim(2,1)*wfx(4)+kim(2,2)*wfy(4)+kim(2,3)*wfz(4);kmvz(4)=kim(3,1)*wfx(4)+kim(3,2)*wfy(4)+kim(3,3)*wfz(4);
                %---- ksi*(mr^-1)*zita* Edge Basis Functions --------------
                kmvwx(1)=ksi(1,1)*mzwx(1)+ksi(1,2)*mzwy(1)+ksi(1,3)*mzwz(1);kmvwy(1)=ksi(2,1)*mzwx(1)+ksi(2,2)*mzwy(1)+ksi(2,3)*mzwz(1);kmvwz(1)=ksi(3,1)*mzwx(1)+ksi(3,2)*mzwy(1)+ksi(3,3)*mzwz(1);
                kmvwx(2)=ksi(1,1)*mzwx(2)+ksi(1,2)*mzwy(2)+ksi(1,3)*mzwz(2);kmvwy(2)=ksi(2,1)*mzwx(2)+ksi(2,2)*mzwy(2)+ksi(2,3)*mzwz(2);kmvwz(2)=ksi(3,1)*mzwx(2)+ksi(3,2)*mzwy(2)+ksi(3,3)*mzwz(2);
                kmvwx(3)=ksi(1,1)*mzwx(3)+ksi(1,2)*mzwy(3)+ksi(1,3)*mzwz(3);kmvwy(3)=ksi(2,1)*mzwx(3)+ksi(2,2)*mzwy(3)+ksi(2,3)*mzwz(3);kmvwz(3)=ksi(3,1)*mzwx(3)+ksi(3,2)*mzwy(3)+ksi(3,3)*mzwz(3);
                kmvwx(4)=ksi(1,1)*mzwx(4)+ksi(1,2)*mzwy(4)+ksi(1,3)*mzwz(4);kmvwy(4)=ksi(2,1)*mzwx(4)+ksi(2,2)*mzwy(4)+ksi(2,3)*mzwz(4);kmvwz(4)=ksi(3,1)*mzwx(4)+ksi(3,2)*mzwy(4)+ksi(3,3)*mzwz(4);
                kmvwx(5)=ksi(1,1)*mzwx(5)+ksi(1,2)*mzwy(5)+ksi(1,3)*mzwz(5);kmvwy(5)=ksi(2,1)*mzwx(5)+ksi(2,2)*mzwy(5)+ksi(2,3)*mzwz(5);kmvwz(5)=ksi(3,1)*mzwx(5)+ksi(3,2)*mzwy(5)+ksi(3,3)*mzwz(5);
                kmvwx(6)=ksi(1,1)*mzwx(6)+ksi(1,2)*mzwy(6)+ksi(1,3)*mzwz(6);kmvwy(6)=ksi(2,1)*mzwx(6)+ksi(2,2)*mzwy(6)+ksi(2,3)*mzwz(6);kmvwz(6)=ksi(3,1)*mzwx(6)+ksi(3,2)*mzwy(6)+ksi(3,3)*mzwz(6);
                %----------------------------------------------------------
                for ii=1:6
                     for jj=1:6,Te(ii,jj)=Te(ii,jj)+Weights(kt)*(wx(ii)*epsilon(1,1)*wx(jj) +wx(ii)*epsilon(1,2)*wy(jj) +wx(ii)*epsilon(1,3)*wz(jj)...
                                                                  +wy(ii)*epsilon(2,1)*wx(jj)+wy(ii)*epsilon(2,2)*wy(jj)+wy(ii)*epsilon(2,3)*wz(jj)...
                                                                  +wz(ii)*epsilon(3,1)*wx(jj)+wz(ii)*epsilon(3,2)*wy(jj)+wz(ii)*epsilon(3,3)*wz(jj));
                                Tam(ii,jj)=Tam(ii,jj)+Weights(kt)*(wx(ii)*kmvwx(jj)+wy(ii)*kmvwy(jj)+wz(ii)*kmvwz(jj));
                                Km(ii,jj)=Km(ii,jj)+Weights(kt)*(wx(ii)*xmwx(jj)+wy(ii)*xmwy(jj)+wz(ii)*xmwz(jj));
                                Pm(ii,jj)=Pm(ii,jj)+Weights(kt)*(rwx(ii)*mzwx(jj)+rwy(ii)*mzwy(jj)+rwz(ii)*mzwz(jj));
                     end
                end
                for ii=1:4,for jj=1:4,Tb(ii,jj)=Tb(ii,jj)+Weights(kt)*(wfx(ii)*wfx(jj)+wfy(ii)*wfy(jj)+wfz(ii)*wfz(jj));end,end
                for ii=1:4
                    for jj=1:6,Fm(ii,jj)=Fm(ii,jj)+Weights(kt)*(wfx(ii)*rwx(jj)+wfy(ii)*rwy(jj)+wfz(ii)*rwz(jj));
                               Fw(ii,jj)=Fw(ii,jj)+Weights(kt)*(wfx(ii)*xwx(jj)+wfy(ii)*xwy(jj)+wfz(ii)*xwz(jj));
                    end
                end
                for ii=1:6
                    for jj=1:4,Am(ii,jj)=Am(ii,jj)+Weights(kt)*(rwx(ii)*mwfx(jj)+rwy(ii)*mwfy(jj)+rwz(ii)*mwfz(jj));
                               Aw(ii,jj)=Aw(ii,jj)+Weights(kt)*(wx(ii)*xwfx(jj)+wy(ii)*xwfy(jj)+wz(ii)*xwfz(jj));
                               Tcm(ii,jj)=Tcm(ii,jj)+Weights(kt)*(wx(ii)*kmvx(jj)+wy(ii)*kmvy(jj)+wz(ii)*kmvz(jj));
                    end
                end
        end,Te=Te*Ve;Tb=Tb*Ve;Am=Am*Ve;Fm=Fm*Ve;Fw=Fw*Ve;Aw=Aw*Ve;Tam=Tam*Ve;Km=Km*Ve;Pm=Pm*Ve;Tcm=Tcm*Ve;
        %------------------------------------------------------------------
        for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                if(boundary.Type=="ABC"),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                   switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end 
                elseif(boundary.Type=="GRA"),Tg = CalculateIntegral(TModel,element,facets(ii),b,c,d,ii);
                    if(boundary.Dispersive),cond=boundary.Param(FreqIndex);else,cond=boundary.Param;end,Tg=Tg.*cond;
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                elseif(boundary.Type=="IBC"),Tbc = CalculateIntegral(TModel,element,facets(ii),b,c,d,ii);
                    if(boundary.Dispersive),cond=boundary.Param(FreqIndex);else,cond=boundary.Param;end,Tbc=Tbc.*cond;
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                end
            end
        end
        %------------------------------------------------------------------
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TEE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TA(counterA)=dsi*dsj*si*sj*Tam(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;P(counterA)=dsi*dsj*si*sj*Pm(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                        
                end
            end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TC(counterA)=si*dsi*sj*dsj*Tcm(ii,jj);
                                   
                end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBB(counterA)=si*sj*dsi*dsj*Tb(ii,jj);end
            end
        end
    end
end
