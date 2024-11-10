%==========================================================================
%{
                E-B Excitation Formulation 
%}
%==========================================================================
function [] = Excitation_Assembly(TModel),AssembledSystem=TModel.Assembled;Frequency=TModel.Frequency;
%----------------------------- Single Frequency Area ----------------------
    if(Frequency.NF==1)
    %------------------------- Dirichlet Boundary Excitation --------------
        if(AssembledSystem.IsDir)
            if(Assembled.Bian),I=zeros(400*TModel.NumberOfElements,1);IA=I;JA=I;TE=I;TB=I;AA=I;FF=I;TS=I;TG=I;TBC=I;P=I;TC=I;TA=I;TP=I;counterA=0;
                               I=zeros(400*TModel.NumberOfElements,1);IB=I;JB=I;B=I;counterB=0;
            else,I=zeros(400*TModel.NumberOfElements,1);IA=I;JA=I;TE=I;TB=I;AA=I;FF=I;TS=I;TG=I;TBC=I;counterA=0;
                 I=zeros(400*TModel.NumberOfElements,1);IB=I;JB=I;B=I;counterB=0;

            end
            %-------------------------------- Domain Assemblies -----------
            for id=1:numel(TModel.Domains),domain=TModel.Domains(id);medium=domain.Medium;
                 switch medium.Type
                     case "Iso" ,[IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,IB,JB,B,counterA,counterB] = Excitation_Dir_Iso(IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,IB,JB,B,counterA,counterB,TModel,domain);
                       case "Anis",[IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,IB,JB,B,counterA,counterB] = Excitation_Dir_Anis(IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,IB,JB,B,counterA,counterB,TModel,domain);
                       case "Bian",[IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,P,TA,TC,IB,JB,B,counterA,counterB] = Excitation_Dir_Bian(IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,P,TA,TC,IB,JB,B,counterA,counterB,TModel,domain);
                  end,disp("domain " + num2str(domain.Index) +"out of " +num2str(numel(TModel.Domains)) + " done");
            end
            %--------------------------------- Sparse Matrices Construction--------
            NonZeros=nnz(IA);IA=IA(1:NonZeros);JA=JA(1:NonZeros);TE=TE(1:NonZeros);TB=TB(1:NonZeros);AA=AA(1:NonZeros);FF=FF(1:NonZeros);TS=TS(1:NonZeros);TG=TG(1:NonZeros);TBC=TBC(1:NonZeros);TP=TP(1:NonZeros);
            NonZerosB=nnz(IB);IB=IB(1:NonZerosB);JB=JB(1:NonZerosB);B=B(1:NonZerosB);
            TE=sparse(IA,JA,TE);TB=sparse(IA,JA,TB);A=sparse(IA,JA,AA);F=sparse(IA,JA,FF);TS=sparse(IA,JA,TS);TG=sparse(IA,JA,TG);TBC=sparse(IA,JA,TBC);TP=soarse(IA,JA,TP);
            B=sparse(IB,JB,B);freq=Frequency.Frequency;
            if(bian_flag),P=P(1:NonZeros);TA=TA(1:NonZeros);TC=TC(1:NonZeros);P=sparse(IA,JA,P);TA=sparse(IA,JA,TA);TC=sparse(IA,JA,TC);
                TModel.Assembled=TModel.Assembled.Excitation(TE,TB,A,F,TS,TG,TBC,TP,P,TA,TC,B,freq);
            else,TModel.Assembled=TModel.Assembled.Excitation(TE,TB,A,F,TS,TG,TBC,TP,freq);
            end
    %-------------------------------- Port Boundary Excitation ------------
        else
            if(Assembled.Bian),I=zeros(400*TModel.NumberOfElements,1);IA=I;JA=I;TEE=I;TBB=I;AA=I;FF=I;TS=I;TG=I;TBC=I;P=I;TC=I;TA=I;TP=I;counterA=0;
            else,I=zeros(400*TModel.NumberOfElements,1);IA=I;JA=I;TEE=I;TBB=I;AA=I;FF=I;TS=I;TG=I;TBC=I;counterA=0;
            end
            %-------------------------------- Domain Assemblies -----------
            for id=1:numel(TModel.Domains),domain=TModel.Domains(id);medium=domain.Medium;
                 switch medium.Type
                       case "Iso" ,[IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,counterA] = Excitation_Port_Iso(IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,counterA,TModel,domain);
                       case "Anis",[IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,counterA] = Excitation_Port_Anis(IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,counterA,TModel,domain);
                       case "Bian",[IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,P,TA,TC,counterA] = Excitation_Port_Bian(IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,P,TA,TC,counterA,TModel,domain);
                  end,disp("domain " + num2str(domain.Index) +"out of " +num2str(numel(TModel.Domains)) + " done");
            end
        end
%------------------------------ Multi Frequency Area ----------------------
    else
    end
end
%======================== Dirichlet Excitation ============================
function [IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,IB,JB,B,counterA,counterB] = Excitation_Dir_Iso(varargin),GaussianQuadratture3D;ElectromagneticConstants;
    if(nargin==17),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};TP=varargin{10};IB=varargin{11};JB=varargin{12};B=varargin{13};
                   counterA=varargin{14};counterB=varargin{15};TModel=varargin{16};domain=varargin(17);
                   medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;freq=TModel.Frequency.Frequency;
    elseif(nargin==18),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};TP=varargin{10};IB=varargin{11};JB=varargin{12};B=varargin{13};
                   counterA=varargin{14};counterB=varargin{15};TModel=varargin{16};domain=varargin(17);FreqIndex=varargin{18};medium=domain.Medium;freq=TModel.Frequency.Frequency(FreqIndex);
                   if(medium.IsDispersive),epsilon=medium.Epsilon(FreqIndex);mu=medium.Mu(FreqIndex);imu=mu^-1;
                   else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;
                   end
    end,k0=2*pi*freq/c0;
    for ie=1:numel(domain.Elements),element=TModel.Elements(domain.Elements(ie));
        b=element.Bs;c=element.Cs;d=element.Ds;Ve=element.Volume;edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];
        edgeSigns=[element.EdgeSigns];facetSigns=[element.FacetSigns];edgeLengths=[edges.Length];facetSurfaces=[facets.Surface];
        Tp=zeros(6,6);Tbc=zeros(6,6);Ts=zeros(6,6);Tg=zeros(6,6);Te=zeros(6,6);Tb=zeros(4,4);Fm=zeros(4,6);Fw=zeros(4,6);Am=zeros(6,4);Aw=zeros(6,4);
        for kt=1:numel(Weights),zeta(1)=Points(1,kt);zeta(2)=Points(2,kt);zeta(3)=Points(3,kt);zeta(4)=Points(4,kt);
            %---------------------- Edge Basis Functions -------------
                 wx(1)=zeta(1)*b(2)-zeta(2)*b(1); wx(2)=zeta(1)*b(3)-zeta(3)*b(1);wx(3)=zeta(1)*b(4)-zeta(4)*b(1);wx(4)=zeta(2)*b(3)-zeta(3)*b(2);wx(5)=zeta(2)*b(4)-zeta(4)*b(2);wx(6)=zeta(3)*b(4)-zeta(4)*b(3);
                 wy(1)=zeta(1)*c(2)-zeta(2)*c(1); wy(2)=zeta(1)*c(3)-zeta(3)*c(1);wy(3)=zeta(1)*c(4)-zeta(4)*c(1);wy(4)=zeta(2)*c(3)-zeta(3)*c(2);wy(5)=zeta(2)*c(4)-zeta(4)*c(2);wy(6)=zeta(3)*c(4)-zeta(4)*c(3);
                 wz(1)=zeta(1)*d(2)-zeta(2)*d(1); wz(2)=zeta(1)*d(3)-zeta(3)*d(1);wz(3)=zeta(1)*d(4)-zeta(4)*d(1);wz(4)=zeta(2)*d(3)-zeta(3)*d(2);wz(5)=zeta(2)*d(4)-zeta(4)*d(2);wz(6)=zeta(3)*d(4)-zeta(4)*d(3);
                 %------------------ Edge Basis Functions Rotations -------
                 
                 %---------------------- Facet Basis Functions ------------
                 wfx(1)=2*zeta(3)*(c(1)*d(2)-c(2)*d(1))+2*zeta(1)*(c(2)*d(3)-c(3)*d(2))+2*zeta(2)*(c(3)*d(1)-c(1)*d(3));wfy(1)=2*zeta(3)*(d(1)*b(2)-d(2)*b(1))+2*zeta(1)*(d(2)*b(3)-d(3)*b(2))+2*zeta(2)*(d(3)*b(1)-d(1)*b(3));wfz(1)=2*zeta(3)*(b(1)*c(2)-b(2)*c(1))+2*zeta(1)*(b(2)*c(3)-b(3)*c(2))+2*zeta(2)*(b(3)*c(1)-b(1)*c(3));
                 wfx(2)=2*zeta(3)*(c(2)*d(4)-c(4)*d(2))+2*zeta(2)*(c(4)*d(3)-c(3)*d(4))+2*zeta(4)*(c(3)*d(2)-c(2)*d(3));wfy(2)=2*zeta(3)*(d(2)*b(4)-d(4)*b(2))+2*zeta(2)*(d(4)*b(3)-d(3)*b(4))+2*zeta(4)*(d(3)*b(2)-d(2)*b(3));wfz(2)=2*zeta(3)*(b(2)*c(4)-b(4)*c(2))+2*zeta(2)*(b(4)*c(3)-b(3)*c(4))+2*zeta(4)*(b(3)*c(2)-b(2)*c(3));
                 wfx(3)=2*zeta(1)*(c(3)*d(4)-c(4)*d(3))+2*zeta(3)*(c(4)*d(1)-c(1)*d(4))+2*zeta(4)*(c(1)*d(3)-c(3)*d(1));wfy(3)=2*zeta(1)*(d(3)*b(4)-d(4)*b(3))+2*zeta(3)*(d(4)*b(1)-d(1)*b(4))+2*zeta(4)*(d(1)*b(3)-d(3)*b(1));wfz(3)=2*zeta(1)*(b(3)*c(4)-b(4)*c(3))+2*zeta(3)*(b(4)*c(1)-b(1)*c(4))+2*zeta(4)*(b(1)*c(3)-b(3)*c(1));
                 wfx(4)=2*zeta(1)*(c(4)*d(2)-c(2)*d(4))+2*zeta(4)*(c(2)*d(1)-c(1)*d(2))+2*zeta(2)*(c(1)*d(4)-c(4)*d(1));wfy(4)=2*zeta(1)*(d(4)*b(2)-d(2)*b(4))+2*zeta(4)*(d(2)*b(1)-d(1)*b(2))+2*zeta(2)*(d(1)*b(4)-d(4)*b(1));wfz(4)=2*zeta(1)*(b(4)*c(2)-b(2)*c(4))+2*zeta(4)*(b(2)*c(1)-b(1)*c(2))+2*zeta(2)*(b(1)*c(4)-b(4)*c(1));
                 %===================== Basis Function Scaling ============
                 switch TModel.AssembledSystem.E_Scaling
                     case 1,wx=wx.*edgeLengths;wy=wy.*edgeLengths;wz=wz.*edgeLengths;rwx=rwx.*edgeLengths;rwy=rwy.*edgeLengths;rwy=rwy.*edgeLengths;rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch TModel.AssembledSystem.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);end
                 end
                                  %==================== Local Matrices =====================
                 for ii=1:6,for jj=1:6,Te(ii,jj)=Te(ii,jj)+Weights(kt)*epsilon*(wx(ii)*wx(jj)+wy(ii)*wy(jj)+wz(ii)*wz(jj));end,end
                 for ii=1:4,for jj=1:4,Tb(ii,jj)=Tb(ii,jj)+Weights(kt)*(wfx(ii)*wfx(jj)+wfy(ii)*wfy(jj)+wfz(ii)*wfz(jj));end,end
                 for ii=1:4,for jj=1:6,Fm(ii,jj)=Fm(ii,jj)+Weights(kt)*(wfx(ii)*rwx(jj)+wfy(ii)*rwy(jj)+wfz(ii)*rwz(jj));end,end
                 for ii=1:6,for jj=1:4,Am(ii,jj)=Am(ii,jj)+Weights(kt)*imu*(rwx(ii)*wfx(jj)+rwy(ii)*wfy(jj)+rwz(ii)*wfz(jj));end,end
        end,Te=Te*Ve;Tb=Tb*Ve;Am=Am*Ve;Fm=Fm*Ve;
     for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                if(boundary.Type=="ABC"),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                    if(medium.IsDispersive),ZW=medium.WaveImpedance(FreqIndex);else,ZW=medium.WaveImpedance;end
                    Ts=1i*c0*m0*(ZW^-1)*Ts;
                   switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                   end 
                elseif(boundary.Type=="ABB"),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                    if(boundary.Dispersive),beta=boundary.Param(FreqIndex);else,beta=boundary.Param;end,k0=2*pi*freq/c0;Ts=Ts*beta*1i/k0;
                     switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                elseif(boundary.Type=="ABZ")
                    if(boundary.Tensor)
                        if(boundary.Dispersive),Z=boundary.Param{FreqIndex};
                        else,Z=boundary.Param;
                        end,Ts=CalculateTensorIntegral(TModel,Z^-1,element,facets(kk),b,c,d,kk);
                    else
                        if(boundary.Dispersive),Z=boundary.Param(FreqIndex);
                        else,Z=boundary.Param;
                        end,Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=Ts*Z^-1;
                    end,Ts=1i*c0*m0*Ts;
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                elseif(boundary.Type=="POR")
                    if(boundary.Dispersive)
                        if(boundary.PortType==0),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param(FreqIndex)/k0)*Ts;Ts=(1i*boundary.Param(FreqIndex)/k0)*Ts;
                        elseif(boundary.PortType==1),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param(FreqIndex)*c0*m0)*Ts;Ts=(1i*boundary.Param(FreqIndex)*c0*m0)*Ts;
                        elseif(boundary.PortType==2),Ts=CalculateTensorIntegral(TModel,boundary.Param{FreqIndex},element,facets(kk),b,c,d,kk);Tp=(2*1i*c0*m0)*Ts;Ts=(1i*c0*m0)*Ts;
                        end
                    else
                        if(boundary.PortType==0),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param/k0)*Ts;Ts=(1i*boundary.Param/k0)*Ts;
                        elseif(boundary.PortType==1),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param*c0*m0)*Ts;Ts=(1i*boundary.Param*c0*m0)*Ts;
                        elseif(boundary.PortType==2),Ts=CalculateTensorIntegral(TModel,boundary.Param,element,facets(kk),b,c,d,kk);Tp=(2*1i*c0*m0)*Ts;Ts=(1i*c0*m0)*Ts;
                        end
                    end
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
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
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TP(counterA)=dsi*dsj*si*sj*Tp(ii,jj);
                end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*k0*Te(ii,jj);end
            end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*Am(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*Am(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*k0*Tb(ii,jj);end
            end
        end  
    end
end
function [IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,IB,JB,B,counterA,counterB] = Excitation_Dir_Anis(varargin),GaussianQuadratture3D;ElectromagneticConstants;
    if(nargin==17),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};TP=varargin{10};IB=varargin{11};JB=varargin{12};B=varargin{13};
                   counterA=varargin{14};counterB=varargin{15};TModel=varargin{16};domain=varargin(17);
                   medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;freq=TModel.Frequency.Frequency;
    elseif(nargin==18),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};TP=varargin{10};IB=varargin{11};JB=varargin{12};B=varargin{13};
                   counterA=varargin{14};counterB=varargin{15};TModel=varargin{16};domain=varargin(17);FreqIndex=varargin{18};medium=domain.Medium;freq=TModel.Frequency.Frequency(FreqIndex);
                   if(medium.IsDispersive),epsilon=medium.Epsilon{FreqIndex};mu=medium.Mu{FreqIndex};imu=mu^-1;
                   else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;
                   end
    end,k0=2*pi*freq/c0;
    for ie=1:numel(domain.Elements),element=TModel.Elements(domain.Elements(ie));disp("element " + num2str(ie) +" out of " + num2str(numel(domain.Elements)));
        b=element.Bs;c=element.Cs;d=element.Ds;Ve=element.Volume;edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];
        edgeSigns=[element.EdgeSigns];facetSigns=[element.FacetSigns];edgeLengths=[edges.Length];facetSurfaces=[facets.Surface];
        Tp=zeros(6,6);Tbc=zeros(6,6);Ts=zeros(6,6);Tg=zeros(6,6);Te=zeros(6,6);Tb=zeros(4,4);Fm=zeros(4,6);Fw=zeros(4,6);Am=zeros(6,4);Aw=zeros(6,4);
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
                 switch TModel.AssembledSystem.E_Scaling
                     case 1,wx=wx.*edgeLengths;wy=wy.*edgeLengths;wz=wz.*edgeLengths;rwx=rwx.*edgeLengths;rwy=rwy.*edgeLengths;rwy=rwy.*edgeLengths;rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch TModel.AssembledSystem.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);end
                 end
                 %------------------------ (mr^-1)*wf----------------------
                 mwfx(1)=imu(1,1)*wfx(1)+imu(1,2)*wfy(1)+imu(1,3)*wfz(1);mwfy(1)=imu(2,1)*wfx(1)+imu(2,2)*wfy(1)+imu(2,3)*wfz(1);mwfz(1)=imu(3,1)*wfx(1)+imu(3,2)*wfy(1)+imu(3,3)*wfz(1);
                 mwfx(2)=imu(1,1)*wfx(2)+imu(1,2)*wfy(2)+imu(1,3)*wfz(2);mwfy(2)=imu(2,1)*wfx(2)+imu(2,2)*wfy(2)+imu(2,3)*wfz(2);mwfz(2)=imu(3,1)*wfx(2)+imu(3,2)*wfy(2)+imu(3,3)*wfz(2);
                 mwfx(3)=imu(1,1)*wfx(3)+imu(1,2)*wfy(3)+imu(1,3)*wfz(3);mwfy(3)=imu(2,1)*wfx(3)+imu(2,2)*wfy(3)+imu(2,3)*wfz(3);mwfz(3)=imu(3,1)*wfx(3)+imu(3,2)*wfy(3)+imu(3,3)*wfz(3);
                 mwfx(4)=imu(1,1)*wfx(4)+imu(1,2)*wfy(4)+imu(1,3)*wfz(4);mwfy(4)=imu(2,1)*wfx(4)+imu(2,2)*wfy(4)+imu(2,3)*wfz(4);mwfz(4)=imu(3,1)*wfx(4)+imu(3,2)*wfy(4)+imu(3,3)*wfz(4);
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
        end,Te=Te*Ve;Tb=Tb*Ve;Am=Am*Ve;Fm=Fm*Ve;
     for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                if(boundary.Type=="ABC"),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                    if(medium.IsDispersive),ZW=medium.WaveImpedance(FreqIndex);else,ZW=medium.WaveImpedance;end
                    Ts=1i*c0*m0*(ZW^-1)*Ts;
                   switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                   end 
                elseif(boundary.Type=="ABB"),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                    if(boundary.Dispersive),beta=boundary.Param(FreqIndex);else,beta=boundary.Param;end,k0=2*pi*freq/c0;Ts=Ts*beta*1i/k0;
                     switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                elseif(boundary.Type=="ABZ")
                    if(boundary.Tensor)
                        if(boundary.Dispersive),Z=boundary.Param{FreqIndex};
                        else,Z=boundary.Param;
                        end,Ts=CalculateTensorIntegral(TModel,Z^-1,element,facets(kk),b,c,d,kk);
                    else
                        if(boundary.Dispersive),Z=boundary.Param(FreqIndex);
                        else,Z=boundary.Param;
                        end,Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=Ts*Z^-1;
                    end,Ts=1i*c0*m0*Ts;
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                elseif(boundary.Type=="POR")
                    if(boundary.Dispersive)
                        if(boundary.PortType==0),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param(FreqIndex)/k0)*Ts;Ts=(1i*boundary.Param(FreqIndex)/k0)*Ts;
                        elseif(boundary.PortType==1),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param(FreqIndex)*c0*m0)*Ts;Ts=(1i*boundary.Param(FreqIndex)*c0*m0)*Ts;
                        elseif(boundary.PortType==2),Ts=CalculateTensorIntegral(TModel,boundary.Param{FreqIndex},element,facets(kk),b,c,d,kk);Tp=(2*1i*c0*m0)*Ts;Ts=(1i*c0*m0)*Ts;
                        end
                    else
                        if(boundary.PortType==0),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param/k0)*Ts;Ts=(1i*boundary.Param/k0)*Ts;
                        elseif(boundary.PortType==1),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param*c0*m0)*Ts;Ts=(1i*boundary.Param*c0*m0)*Ts;
                        elseif(boundary.PortType==2),Ts=CalculateTensorIntegral(TModel,boundary.Param,element,facets(kk),b,c,d,kk);Tp=(2*1i*c0*m0)*Ts;Ts=(1i*c0*m0)*Ts;
                        end
                    end
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
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
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TP(counterA)=dsi*dsj*si*sj*Tp(ii,jj);
                end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*k0*Te(ii,jj);end
            end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*Am(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*Am(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*k0*Tb(ii,jj);end
            end
        end  
    end
end
function [IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,P,TA,TC,IB,JB,B,counterA,counterB] = Excitation_Dir_Bian(varargin),GaussianQuadratture3D;
    if(nargin==20),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};TP=varargin{10};P=varargin{11};TA=varargin{12};TC=varargin{13};IB=varargin{14};
                   JB=varargin{15};B=varargin{16};counterA=varargin{17};
                   counterB=varargin{18};TModel=varargin{19};domain=varargin{20};
                   medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;ksi=medium.Ksi;zita=medium.Zita;
    elseif(nargin==21),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};TP=varargin{10};P=varargin{11};TA=varargin{12};TC=varargin{13};IB=varargin{14};
                   JB=varargin{15};B=varargin{16};counterA=varargin{17};counterB=varargin{18};TModel=varargin{19};domain=varargin{20};
                       FreqIndex=varargin{21};medium=domain.Medium;
                       if(medium.IsDispersive),epsilon=medium.Epsilon{FreqIndex};mu=medium.Mu{FreqIndex};imu=mu^-1;ksi=medium.Ksi{FreqIndex};zita=medium.Zita{FreqIndex};
                       else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;ksi=medium.Ksi;zita=medium.Zita;
                       end
    end,imz=imu*zita;kim=ksi*imu;
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
                 switch TModel.AssembledSystem.E_Scaling
                     case 1,wx=wx.*edgeLengths;wy=wy.*edgeLengths;wz=wz.*edgeLengths;rwx=rwx.*edgeLengths;rwy=rwy.*edgeLengths;rwy=rwy.*edgeLengths;rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch TModel.AssembledSystem.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);end
                 end
                 %-------------  (mr^-1) * Facet Basis Functions-----------
                 mwfx(1)=imu(1,1)*wfx(1)+imu(1,2)*wfy(1)+imu(1,3)*wfz(1);mwfy(1)=imu(2,1)*wfx(1)+imu(2,2)*wfy(1)+imu(2,3)*wfz(1);mwfz(1)=imu(3,1)*wfx(1)+imu(3,2)*wfy(1)+imu(3,3)*wfz(1);
                 mwfx(2)=imu(1,1)*wfx(2)+imu(1,2)*wfy(2)+imu(1,3)*wfz(2);mwfy(2)=imu(2,1)*wfx(2)+imu(2,2)*wfy(2)+imu(2,3)*wfz(2);mwfz(2)=imu(3,1)*wfx(2)+imu(3,2)*wfy(2)+imu(3,3)*wfz(2);
                 mwfx(3)=imu(1,1)*wfx(3)+imu(1,2)*wfy(3)+imu(1,3)*wfz(3);mwfy(3)=imu(2,1)*wfx(3)+imu(2,2)*wfy(3)+imu(2,3)*wfz(3);mwfz(3)=imu(3,1)*wfx(3)+imu(3,2)*wfy(3)+imu(3,3)*wfz(3);
                 mwfx(4)=imu(1,1)*wfx(4)+imu(1,2)*wfy(4)+imu(1,3)*wfz(4);mwfy(4)=imu(2,1)*wfx(4)+imu(2,2)*wfy(4)+imu(2,3)*wfz(4);mwfz(4)=imu(3,1)*wfx(4)+imu(3,2)*wfy(4)+imu(3,3)*wfz(4);
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
               %-------------- x cross Edge Basis Functions -------------
                xwx(1)=0;xwy(1)=-wz(1);xwz(1)=wy(1);
                xwx(2)=0;xwy(2)=-wz(2);xwz(2)=wy(2);
                xwx(3)=0;xwy(3)=-wz(3);xwz(3)=wy(3);
                xwx(4)=0;xwy(4)=-wz(4);xwz(4)=wy(4);
                xwx(5)=0;xwy(5)=-wz(5);xwz(5)=wy(5);
                xwx(6)=0;xwy(6)=-wz(6);xwz(6)=wy(6);
                %-------------- x cross mr^-1 Facet Basis Functions ------
                xwfx(1)=0;xwfy(1)=-mwfz(1);xwfz(1)=mwfy(1);
                xwfx(2)=0;xwfy(2)=-mwfz(2);xwfz(2)=mwfy(2);
                xwfx(3)=0;xwfy(3)=-mwfz(3);xwfz(3)=mwfy(3);
                xwfx(4)=0;xwfy(4)=-mwfz(4);xwfz(4)=mwfy(4);
                %-------------- x cross (mr^-1)*zita* Edge Basis Functions 
                xmwx(1)=0;xmwy(1)=-mzwz(1);xmwz(1)=mzwy(1);
                xmwx(2)=0;xmwy(2)=-mzwz(2);xmwz(2)=mzwy(2);
                xmwx(3)=0;xmwy(3)=-mzwz(3);xmwz(3)=mzwy(3);
                xmwx(4)=0;xmwy(4)=-mzwz(4);xmwz(4)=mzwy(4);
                xmwx(5)=0;xmwy(5)=-mzwz(5);xmwz(5)=mzwy(5);
                xmwx(6)=0;xmwy(6)=-mzwz(6);xmwz(6)=mzwy(6);
                %----------------------------------------------------------
                for ii=1:6
                     for jj=1:6,Te(ii,jj)=Te(ii,jj)+Weights(kt)*(wx(ii)*epsilon(1,1)*wx(jj) +wx(ii)*epsilon(1,2)*wy(jj) +wx(ii)*epsilon(1,3)*wz(jj)...
                                                                  +wy(ii)*epsilon(2,1)*wx(jj)+wy(ii)*epsilon(2,2)*wy(jj)+wy(ii)*epsilon(2,3)*wz(jj)...
                                                                  +wz(ii)*epsilon(3,1)*wx(jj)+wz(ii)*epsilon(3,2)*wy(jj)+wz(ii)*epsilon(3,3)*wz(jj));
                                Tam(ii,jj)=Tam(ii,jj)+Weights(kt)*(wx(ii)*kmvwx(jj)+wy(ii)*kmvwy(jj)+wz(ii)*kmvwz(jj));
                                Pm(ii,jj)=Pm(ii,jj)+Weights(kt)*(rwx(ii)*mzwx(jj)+rwy(ii)*mzwy(jj)+rwz(ii)*mzwz(jj));
                     end
                end
                for ii=1:4,for jj=1:4,Tb(ii,jj)=Tb(ii,jj)+Weights(kt)*(wfx(ii)*wfx(jj)+wfy(ii)*wfy(jj)+wfz(ii)*wfz(jj));end,end
                for ii=1:4,for jj=1:6,Fm(ii,jj)=Fm(ii,jj)+Weights(kt)*(wfx(ii)*rwx(jj)+wfy(ii)*rwy(jj)+wfz(ii)*rwz(jj));end,end
                for ii=1:6
                    for jj=1:4,Am(ii,jj)=Am(ii,jj)+Weights(kt)*(rwx(ii)*mwfx(jj)+rwy(ii)*mwfy(jj)+rwz(ii)*mwfz(jj));
                               Tcm(ii,jj)=Tcm(ii,jj)+Weights(kt)*(wx(ii)*kmvx(jj)+wy(ii)*kmvy(jj)+wz(ii)*kmvz(jj));
                    end
                end
        end,Te=Te*Ve;Tb=Tb*Ve;Am=Am*Ve;Fm=Fm*Ve;Tam=Tam*Ve;Pm=Pm*Ve;Tcm=Tcm*Ve;
        %------------------------------------------------------------------
        for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                if(boundary.Type=="ABC")
                    if(medium.IsDispersive),ZW=medium.WaveImpedance{FreqIndex};else,ZW=medium.WaveImpedance;end
                    Ts = CalculateTensorIntegral(TModel,ZW^-1,element,facets(kk),b,c,d,kk);
                    Ts=1i*c0*m0*Ts;
                   switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                   end 
                elseif(boundary.Type=="ABB"),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                    if(boundary.Dispersive),beta=boundary.Param(FreqIndex);else,beta=boundary.Param;end,Ts=Ts*beta*1i/(2*pi*freq/c0);
                     switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                elseif(boundary.Type=="ABZ")
                    if(boundary.Tensor)
                        if(boundary.Dispersive),Z=boundary.Param{FreqIndex};
                        else,Z=boundary.Param;
                        end,Ts=CalculateTensorIntegral(TModel,Z^-1,element,facets(kk),b,c,d,kk);
                    else
                        if(boundary.Dispersive),Z=boundary.Param(FreqIndex);
                        else,Z=boundary.Param;
                        end,Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=Ts*Z^-1;
                    end,Ts=1i*c0*m0*Ts;
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end 
                elseif(boundary.Type=="POR")
                    if(boundary.Dispersive)
                        if(boundary.PortType==0),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param(FreqIndex)/k0)*Ts;Ts=(1i*boundary.Param(FreqIndex)/k0)*Ts;
                        elseif(boundary.PortType==1),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param(FreqIndex)*c0*m0)*Ts;Ts=(1i*boundary.Param(FreqIndex)*c0*m0)*Ts;
                        elseif(boundary.PortType==2),Ts=CalculateTensorIntegral(TModel,boundary.Param{FreqIndex},element,facets(kk),b,c,d,kk);Tp=(2*1i*c0*m0)*Ts;Ts=(1i*c0*m0)*Ts;
                        end
                    else
                        if(boundary.PortType==0),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param/k0)*Ts;Ts=(1i*boundary.Param/k0)*Ts;
                        elseif(boundary.PortType==1),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param*c0*m0)*Ts;Ts=(1i*boundary.Param*c0*m0)*Ts;
                        elseif(boundary.PortType==2),Ts=CalculateTensorIntegral(TModel,boundary.Param,element,facets(kk),b,c,d,kk);Tp=(2*1i*c0*m0)*Ts;Ts=(1i*c0*m0)*Ts;
                        end
                    end
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
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
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TA(counterA)=dsi*dsj*si*sj*Tam(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;P(counterA)=dsi*dsj*si*sj*Pm(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=dsi*dsj*si*sj*k0*Te(ii,jj);
                                   counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=-dsi*dsj*si*sj*2*pi*freq*c0*Tam(ii,jj);
                                   counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=-1i*c0*dsi*dsj*si*sj*Pm(ii,jj);
                end
            end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TC(counterA)=si*dsi*sj*dsj*Tcm(ii,jj);
                end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterA)=kj;B(counterB)=si*dsi*sj*dsj*Am(ii,jj);
                                   counterB=counterB+1;IB(counterB)=di;JB(counterA)=kj;B(counterB)=-1i*2*pi*freq*si*dsi*sj*dsj*Tcm(ii,jj);

                end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;B(counterB)=si*dsi*sj*dsj*Fm(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IA(counterB)=di;JB(counterB)=kj;B(counterB)=si*sj*dsi*dsj*k0*Tb(ii,jj);end
            end
        end
    end
end
%======================== Port Boundary Excitataion =======================
function [IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,counterA] = Excitation_Port_Iso(varargin),GaussianQuadratture3D;ElectromagneticConstants;
    if(nargin==13),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};TP=varargin{10};
                   counterA=varargin{11};TModel=varargin{12};domain=varargin(13);
                   medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;freq=TModel.Frequency.Frequency;
    elseif(nargin==14),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};TP=varargin{10};
                   counterA=varargin{11};TModel=varargin{12};domain=varargin(13);FreqIndex=varargin{14};medium=domain.Medium;freq=TModel.Frequency.Frequency(FreqIndex);
                   if(medium.IsDispersive),epsilon=medium.Epsilon(FreqIndex);mu=medium.Mu(FreqIndex);imu=mu^-1;
                   else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;
                   end
    end,k0=2*pi*freq/c0;
    for ie=1:numel(domain.Elements),element=TModel.Elements(domain.Elements(ie));
        b=element.Bs;c=element.Cs;d=element.Ds;Ve=element.Volume;edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];
        edgeSigns=[element.EdgeSigns];facetSigns=[element.FacetSigns];edgeLengths=[edges.Length];facetSurfaces=[facets.Surface];
        Tp=zeros(6,6);Tbc=zeros(6,6);Ts=zeros(6,6);Tg=zeros(6,6);Te=zeros(6,6);Tb=zeros(4,4);Fm=zeros(4,6);Fw=zeros(4,6);Am=zeros(6,4);Aw=zeros(6,4);
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
                 switch TModel.AssembledSystem.E_Scaling
                     case 1,wx=wx.*edgeLengths;wy=wy.*edgeLengths;wz=wz.*edgeLengths;rwx=rwx.*edgeLengths;rwy=rwy.*edgeLengths;rwy=rwy.*edgeLengths;rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch TModel.AssembledSystem.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);end
                 end
                                  %==================== Local Matrices =====================
                 for ii=1:6,for jj=1:6,Te(ii,jj)=Te(ii,jj)+Weights(kt)*epsilon*(wx(ii)*wx(jj)+wy(ii)*wy(jj)+wz(ii)*wz(jj));end,end
                 for ii=1:4,for jj=1:4,Tb(ii,jj)=Tb(ii,jj)+Weights(kt)*(wfx(ii)*wfx(jj)+wfy(ii)*wfy(jj)+wfz(ii)*wfz(jj));end,end
                 for ii=1:4,for jj=1:6,Fm(ii,jj)=Fm(ii,jj)+Weights(kt)*(wfx(ii)*rwx(jj)+wfy(ii)*rwy(jj)+wfz(ii)*rwz(jj));end,end
                 for ii=1:6,for jj=1:4,Am(ii,jj)=Am(ii,jj)+Weights(kt)*imu*(rwx(ii)*wfx(jj)+rwy(ii)*wfy(jj)+rwz(ii)*wfz(jj));end,end
        end,Te=Te*Ve;Tb=Tb*Ve;Am=Am*Ve;Fm=Fm*Ve;
     for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                if(boundary.Type=="ABC"),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                    if(medium.IsDispersive),ZW=medium.WaveImpedance(FreqIndex);else,ZW=medium.WaveImpedance;end
                    Ts=1i*c0*m0*(ZW^-1)*Ts;
                   switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                   end 
                elseif(boundary.Type=="ABB"),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                    if(boundary.Dispersive),beta=boundary.Param(FreqIndex);else,beta=boundary.Param;end,k0=2*pi*freq/c0;Ts=Ts*beta*1i/k0;
                     switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                elseif(boundary.Type=="ABZ")
                    if(boundary.Tensor)
                        if(boundary.Dispersive),Z=boundary.Param{FreqIndex};
                        else,Z=boundary.Param;
                        end,Ts=CalculateTensorIntegral(TModel,Z^-1,element,facets(kk),b,c,d,kk);
                    else
                        if(boundary.Dispersive),Z=boundary.Param(FreqIndex);
                        else,Z=boundary.Param;
                        end,Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=Ts*Z^-1;
                    end,Ts=1i*c0*m0*Ts;
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                elseif(boundary.Type=="POR")
                    if(boundary.Dispersive)
                        if(boundary.PortType==0),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param(FreqIndex)/k0)*Ts;Ts=(1i*boundary.Param(FreqIndex)/k0)*Ts;
                        elseif(boundary.PortType==1),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param(FreqIndex)*c0*m0)*Ts;Ts=(1i*boundary.Param(FreqIndex)*c0*m0)*Ts;
                        elseif(boundary.PortType==2),Ts=CalculateTensorIntegral(TModel,boundary.Param{FreqIndex},element,facets(kk),b,c,d,kk);Tp=(2*1i*c0*m0)*Ts;Ts=(1i*c0*m0)*Ts;
                        end
                    else
                        if(boundary.PortType==0),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param/k0)*Ts;Ts=(1i*boundary.Param/k0)*Ts;
                        elseif(boundary.PortType==1),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param*c0*m0)*Ts;Ts=(1i*boundary.Param*c0*m0)*Ts;
                        elseif(boundary.PortType==2),Ts=CalculateTensorIntegral(TModel,boundary.Param,element,facets(kk),b,c,d,kk);Tp=(2*1i*c0*m0)*Ts;Ts=(1i*c0*m0)*Ts;
                        end
                    end
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
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
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TP(counterA)=dsi*dsj*si*sj*Tp(ii,jj);
                end
            end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
             end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj);end
            end
        end  
    end
end
function [IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,counterA] = Excitation_Port_Anis(varargin),GaussianQuadratture3D;ElectromagneticConstants;
    if(nargin==13),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};TP=varargin{10};counterA=varargin{11};TModel=varargin{12};domain=varargin(13);
                   medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;freq=TModel.Frequency.Frequency;
    elseif(nargin==14),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};TP=varargin{10};counterA=varargin{11};TModel=varargin{12};domain=varargin(13);FreqIndex=varargin{14};medium=domain.Medium;freq=TModel.Frequency.Frequency(FreqIndex);
                   if(medium.IsDispersive),epsilon=medium.Epsilon{FreqIndex};mu=medium.Mu{FreqIndex};imu=mu^-1;
                   else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;
                   end
    end,k0=2*pi*freq/c0;
    for ie=1:numel(domain.Elements),element=TModel.Elements(domain.Elements(ie));disp("element " + num2str(ie) +" out of " + num2str(numel(domain.Elements)));
        b=element.Bs;c=element.Cs;d=element.Ds;Ve=element.Volume;edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];
        edgeSigns=[element.EdgeSigns];facetSigns=[element.FacetSigns];edgeLengths=[edges.Length];facetSurfaces=[facets.Surface];
        Tp=zeros(6,6);Tbc=zeros(6,6);Ts=zeros(6,6);Tg=zeros(6,6);Te=zeros(6,6);Tb=zeros(4,4);Fm=zeros(4,6);Fw=zeros(4,6);Am=zeros(6,4);Aw=zeros(6,4);
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
                 switch TModel.AssembledSystem.E_Scaling
                     case 1,wx=wx.*edgeLengths;wy=wy.*edgeLengths;wz=wz.*edgeLengths;rwx=rwx.*edgeLengths;rwy=rwy.*edgeLengths;rwy=rwy.*edgeLengths;rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch TModel.AssembledSystem.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);end
                 end
                 %------------------------ (mr^-1)*wf----------------------
                 mwfx(1)=imu(1,1)*wfx(1)+imu(1,2)*wfy(1)+imu(1,3)*wfz(1);mwfy(1)=imu(2,1)*wfx(1)+imu(2,2)*wfy(1)+imu(2,3)*wfz(1);mwfz(1)=imu(3,1)*wfx(1)+imu(3,2)*wfy(1)+imu(3,3)*wfz(1);
                 mwfx(2)=imu(1,1)*wfx(2)+imu(1,2)*wfy(2)+imu(1,3)*wfz(2);mwfy(2)=imu(2,1)*wfx(2)+imu(2,2)*wfy(2)+imu(2,3)*wfz(2);mwfz(2)=imu(3,1)*wfx(2)+imu(3,2)*wfy(2)+imu(3,3)*wfz(2);
                 mwfx(3)=imu(1,1)*wfx(3)+imu(1,2)*wfy(3)+imu(1,3)*wfz(3);mwfy(3)=imu(2,1)*wfx(3)+imu(2,2)*wfy(3)+imu(2,3)*wfz(3);mwfz(3)=imu(3,1)*wfx(3)+imu(3,2)*wfy(3)+imu(3,3)*wfz(3);
                 mwfx(4)=imu(1,1)*wfx(4)+imu(1,2)*wfy(4)+imu(1,3)*wfz(4);mwfy(4)=imu(2,1)*wfx(4)+imu(2,2)*wfy(4)+imu(2,3)*wfz(4);mwfz(4)=imu(3,1)*wfx(4)+imu(3,2)*wfy(4)+imu(3,3)*wfz(4);
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
        end,Te=Te*Ve;Tb=Tb*Ve;Am=Am*Ve;Fm=Fm*Ve;
     for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                if(boundary.Type=="ABC"),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                    if(medium.IsDispersive),ZW=medium.WaveImpedance(FreqIndex);else,ZW=medium.WaveImpedance;end
                    Ts=1i*c0*m0*(ZW^-1)*Ts;
                   switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                   end 
                elseif(boundary.Type=="ABB"),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                    if(boundary.Dispersive),beta=boundary.Param(FreqIndex);else,beta=boundary.Param;end,k0=2*pi*freq/c0;Ts=Ts*beta*1i/k0;
                     switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                elseif(boundary.Type=="ABZ")
                    if(boundary.Tensor)
                        if(boundary.Dispersive),Z=boundary.Param{FreqIndex};
                        else,Z=boundary.Param;
                        end,Ts=CalculateTensorIntegral(TModel,Z^-1,element,facets(kk),b,c,d,kk);
                    else
                        if(boundary.Dispersive),Z=boundary.Param(FreqIndex);
                        else,Z=boundary.Param;
                        end,Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=Ts*Z^-1;
                    end,Ts=1i*c0*m0*Ts;
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                elseif(boundary.Type=="POR")
                    if(boundary.Dispersive)
                        if(boundary.PortType==0),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param(FreqIndex)/k0)*Ts;Ts=(1i*boundary.Param(FreqIndex)/k0)*Ts;
                        elseif(boundary.PortType==1),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param(FreqIndex)*c0*m0)*Ts;Ts=(1i*boundary.Param(FreqIndex)*c0*m0)*Ts;
                        elseif(boundary.PortType==2),Ts=CalculateTensorIntegral(TModel,boundary.Param{FreqIndex},element,facets(kk),b,c,d,kk);Tp=(2*1i*c0*m0)*Ts;Ts=(1i*c0*m0)*Ts;
                        end
                    else
                        if(boundary.PortType==0),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param/k0)*Ts;Ts=(1i*boundary.Param/k0)*Ts;
                        elseif(boundary.PortType==1),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param*c0*m0)*Ts;Ts=(1i*boundary.Param*c0*m0)*Ts;
                        elseif(boundary.PortType==2),Ts=CalculateTensorIntegral(TModel,boundary.Param,element,facets(kk),b,c,d,kk);Tp=(2*1i*c0*m0)*Ts;Ts=(1i*c0*m0)*Ts;
                        end
                    end
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
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
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TP(counterA)=dsi*dsj*si*sj*Tp(ii,jj);
                end
             end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj);end
            end
        end  
    end
end
function [IA,JA,TE,TB,AA,FF,TS,TG,TBC,TP,P,TA,TC,counterA] = Excitation_Port_Bian(varargin),GaussianQuadratture3D;
    if(nargin==16),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};TP=varargin{10};P=varargin{11};TA=varargin{12};TC=varargin{13};counterA=varargin{14};
                   TModel=varargin{15};domain=varargin{16};
                   medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;ksi=medium.Ksi;zita=medium.Zita;
    elseif(nargin==17),IA=varargin{1};JA=varargin{2};TE=varargin{3};TB=varargin{4};AA=varargin{5};FF=varargin{6};
                   TS=varargin{7};TG=varargin{8};TBC=varargin{9};TP=varargin{10};P=varargin{11};TA=varargin{12};TC=varargin{13};counterA=varargin{14};
                   TModel=varargin{15};domain=varargin{16};
                       FreqIndex=varargin{17};medium=domain.Medium;
                       if(medium.IsDispersive),epsilon=medium.Epsilon{FreqIndex};mu=medium.Mu{FreqIndex};imu=mu^-1;ksi=medium.Ksi{FreqIndex};zita=medium.Zita{FreqIndex};
                       else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;ksi=medium.Ksi;zita=medium.Zita;
                       end
    end,imz=imu*zita;kim=ksi*imu;
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
                 switch TModel.AssembledSystem.E_Scaling
                     case 1,wx=wx.*edgeLengths;wy=wy.*edgeLengths;wz=wz.*edgeLengths;rwx=rwx.*edgeLengths;rwy=rwy.*edgeLengths;rwy=rwy.*edgeLengths;rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*AssembledSystem.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch TModel.AssembledSystem.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*AssembledSystem.ScalingVec_B(facets(ii).Index);end
                 end
                 %-------------  (mr^-1) * Facet Basis Functions-----------
                 mwfx(1)=imu(1,1)*wfx(1)+imu(1,2)*wfy(1)+imu(1,3)*wfz(1);mwfy(1)=imu(2,1)*wfx(1)+imu(2,2)*wfy(1)+imu(2,3)*wfz(1);mwfz(1)=imu(3,1)*wfx(1)+imu(3,2)*wfy(1)+imu(3,3)*wfz(1);
                 mwfx(2)=imu(1,1)*wfx(2)+imu(1,2)*wfy(2)+imu(1,3)*wfz(2);mwfy(2)=imu(2,1)*wfx(2)+imu(2,2)*wfy(2)+imu(2,3)*wfz(2);mwfz(2)=imu(3,1)*wfx(2)+imu(3,2)*wfy(2)+imu(3,3)*wfz(2);
                 mwfx(3)=imu(1,1)*wfx(3)+imu(1,2)*wfy(3)+imu(1,3)*wfz(3);mwfy(3)=imu(2,1)*wfx(3)+imu(2,2)*wfy(3)+imu(2,3)*wfz(3);mwfz(3)=imu(3,1)*wfx(3)+imu(3,2)*wfy(3)+imu(3,3)*wfz(3);
                 mwfx(4)=imu(1,1)*wfx(4)+imu(1,2)*wfy(4)+imu(1,3)*wfz(4);mwfy(4)=imu(2,1)*wfx(4)+imu(2,2)*wfy(4)+imu(2,3)*wfz(4);mwfz(4)=imu(3,1)*wfx(4)+imu(3,2)*wfy(4)+imu(3,3)*wfz(4);
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
               %-------------- x cross Edge Basis Functions -------------
                xwx(1)=0;xwy(1)=-wz(1);xwz(1)=wy(1);
                xwx(2)=0;xwy(2)=-wz(2);xwz(2)=wy(2);
                xwx(3)=0;xwy(3)=-wz(3);xwz(3)=wy(3);
                xwx(4)=0;xwy(4)=-wz(4);xwz(4)=wy(4);
                xwx(5)=0;xwy(5)=-wz(5);xwz(5)=wy(5);
                xwx(6)=0;xwy(6)=-wz(6);xwz(6)=wy(6);
                %-------------- x cross mr^-1 Facet Basis Functions ------
                xwfx(1)=0;xwfy(1)=-mwfz(1);xwfz(1)=mwfy(1);
                xwfx(2)=0;xwfy(2)=-mwfz(2);xwfz(2)=mwfy(2);
                xwfx(3)=0;xwfy(3)=-mwfz(3);xwfz(3)=mwfy(3);
                xwfx(4)=0;xwfy(4)=-mwfz(4);xwfz(4)=mwfy(4);
                %-------------- x cross (mr^-1)*zita* Edge Basis Functions 
                xmwx(1)=0;xmwy(1)=-mzwz(1);xmwz(1)=mzwy(1);
                xmwx(2)=0;xmwy(2)=-mzwz(2);xmwz(2)=mzwy(2);
                xmwx(3)=0;xmwy(3)=-mzwz(3);xmwz(3)=mzwy(3);
                xmwx(4)=0;xmwy(4)=-mzwz(4);xmwz(4)=mzwy(4);
                xmwx(5)=0;xmwy(5)=-mzwz(5);xmwz(5)=mzwy(5);
                xmwx(6)=0;xmwy(6)=-mzwz(6);xmwz(6)=mzwy(6);
                %----------------------------------------------------------
                for ii=1:6
                     for jj=1:6,Te(ii,jj)=Te(ii,jj)+Weights(kt)*(wx(ii)*epsilon(1,1)*wx(jj) +wx(ii)*epsilon(1,2)*wy(jj) +wx(ii)*epsilon(1,3)*wz(jj)...
                                                                  +wy(ii)*epsilon(2,1)*wx(jj)+wy(ii)*epsilon(2,2)*wy(jj)+wy(ii)*epsilon(2,3)*wz(jj)...
                                                                  +wz(ii)*epsilon(3,1)*wx(jj)+wz(ii)*epsilon(3,2)*wy(jj)+wz(ii)*epsilon(3,3)*wz(jj));
                                Tam(ii,jj)=Tam(ii,jj)+Weights(kt)*(wx(ii)*kmvwx(jj)+wy(ii)*kmvwy(jj)+wz(ii)*kmvwz(jj));
                                Pm(ii,jj)=Pm(ii,jj)+Weights(kt)*(rwx(ii)*mzwx(jj)+rwy(ii)*mzwy(jj)+rwz(ii)*mzwz(jj));
                     end
                end
                for ii=1:4,for jj=1:4,Tb(ii,jj)=Tb(ii,jj)+Weights(kt)*(wfx(ii)*wfx(jj)+wfy(ii)*wfy(jj)+wfz(ii)*wfz(jj));end,end
                for ii=1:4,for jj=1:6,Fm(ii,jj)=Fm(ii,jj)+Weights(kt)*(wfx(ii)*rwx(jj)+wfy(ii)*rwy(jj)+wfz(ii)*rwz(jj));end,end
                for ii=1:6
                    for jj=1:4,Am(ii,jj)=Am(ii,jj)+Weights(kt)*(rwx(ii)*mwfx(jj)+rwy(ii)*mwfy(jj)+rwz(ii)*mwfz(jj));
                               Tcm(ii,jj)=Tcm(ii,jj)+Weights(kt)*(wx(ii)*kmvx(jj)+wy(ii)*kmvy(jj)+wz(ii)*kmvz(jj));
                    end
                end
        end,Te=Te*Ve;Tb=Tb*Ve;Am=Am*Ve;Fm=Fm*Ve;Tam=Tam*Ve;Pm=Pm*Ve;Tcm=Tcm*Ve;
        %------------------------------------------------------------------
        for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                if(boundary.Type=="ABC")
                    if(medium.IsDispersive),ZW=medium.WaveImpedance{FreqIndex};else,ZW=medium.WaveImpedance;end
                    Ts = CalculateTensorIntegral(TModel,ZW^-1,element,facets(kk),b,c,d,kk);
                    Ts=1i*c0*m0*Ts;
                   switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                   end 
                elseif(boundary.Type=="ABB"),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                    if(boundary.Dispersive),beta=boundary.Param(FreqIndex);else,beta=boundary.Param;end,Ts=Ts*beta*1i/(2*pi*freq/c0);
                     switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end
                elseif(boundary.Type=="ABZ")
                    if(boundary.Tensor)
                        if(boundary.Dispersive),Z=boundary.Param{FreqIndex};
                        else,Z=boundary.Param;
                        end,Ts=CalculateTensorIntegral(TModel,Z^-1,element,facets(kk),b,c,d,kk);
                    else
                        if(boundary.Dispersive),Z=boundary.Param(FreqIndex);
                        else,Z=boundary.Param;
                        end,Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=Ts*Z^-1;
                    end,Ts=1i*c0*m0*Ts;
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
                    end 
                elseif(boundary.Type=="POR")
                    if(boundary.Dispersive)
                        if(boundary.PortType==0),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param(FreqIndex)/k0)*Ts;Ts=(1i*boundary.Param(FreqIndex)/k0)*Ts;
                        elseif(boundary.PortType==1),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param(FreqIndex)*c0*m0)*Ts;Ts=(1i*boundary.Param(FreqIndex)*c0*m0)*Ts;
                        elseif(boundary.PortType==2),Ts=CalculateTensorIntegral(TModel,boundary.Param{FreqIndex},element,facets(kk),b,c,d,kk);Tp=(2*1i*c0*m0)*Ts;Ts=(1i*c0*m0)*Ts;
                        end
                    else
                        if(boundary.PortType==0),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param/k0)*Ts;Ts=(1i*boundary.Param/k0)*Ts;
                        elseif(boundary.PortType==1),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Tp=(2*1i*boundary.Param*c0*m0)*Ts;Ts=(1i*boundary.Param*c0*m0)*Ts;
                        elseif(boundary.PortType==2),Ts=CalculateTensorIntegral(TModel,boundary.Param,element,facets(kk),b,c,d,kk);Tp=(2*1i*c0*m0)*Ts;Ts=(1i*c0*m0)*Ts;
                        end
                    end
                    switch AssembledSystem.E_Scaling
                        case 1,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                        case 2,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*AssembledSystem.ScalingVec_E(edges(ii).Index)*AssembledSystem.ScalingVec_E(edges(jj).Index);end,end
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
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TA(counterA)=dsi*dsj*si*sj*Tam(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;P(counterA)=dsi*dsj*si*sj*Pm(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                end
            end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(ii).KnownIndex;
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
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(ii).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj);end
            end
        end
    end
end
%--------------------------------------------------------------------------
function [SurfaceIntegral] = CalculateIntegral(TModel,element,facet,b,c,d,facetIndex),nvec=findNormaVector(TModel,element,facetIndex);SurfaceIntegral=zeros(6,6);Ae=facet.Surface;
    switch facetIndex
        case 1,dz1=[b(1);c(1);d(1)];dz2=[b(2);c(2);d(2)];dz3=[b(3);c(3);d(3)];
               dt1=cross(nvec,dz1);dt1=-cross(nvec,dt1);dt2=cross(nvec,dz2);dt2=-cross(nvec,dt2);dt3=cross(nvec,dz3);dt3=-cross(nvec,dt3);
               SurfaceIntegral(1,1)=Ae*(dt2'*dt2 + dt1'*dt1)/6 - Ae*(dt2'*dt1 + dt1'*dt2)/12;
               SurfaceIntegral(1,2)=Ae*(dt2'*dt3)/6+Ae*(dt1'*dt1 - dt2'*dt1 - dt1'*dt2)/12;
               SurfaceIntegral(1,4)=-Ae*(dt1'*dt3)/6 +Ae*(dt1'*dt1 - dt2'*dt1 - dt1'*dt3)/12;
               SurfaceIntegral(2,2)=Ae*(dt3'*dt3 + dt1'*dt1)/6 -Ae*(dt3'*dt1 + dt1'*dt3)/12;
               SurfaceIntegral(2,1)=Ae*(dt3'*dt2)/6 + Ae*(dt1'*dt1 -dt3'*dt1 - dt1'*dt2)/12;
               SurfaceIntegral(2,4)=Ae*(dt1'*dt2)/6 + Ae*(dt3'*dt3 - dt3'*dt2 -dt1'*dt3)/12;
               SurfaceIntegral(4,4)=Ae*(dt3'*dt3 +dt2'*dt2)/6 - Ae*(dt3'*dt2+ dt2'*dt3)/12;
               SurfaceIntegral(4,1)=-Ae*(dt3'*dt1)/6 + Ae*(dt3'*dt2 +dt2'*dt1 -dt2'*dt2)/12;
               SurfaceIntegral(4,2)=Ae*(dt2'*dt1)/6 +Ae*(dt3'*dt3 -dt3'*dt1 -dt2'*dt3)/12;
        case 2,dz2=[b(2);c(2);d(2)];dz4=[b(4);c(4);d(4)];dz3=[b(3);c(3);d(3)];
               dt2=cross(nvec,dz2);dt2=-cross(nvec,dt2);dt4=cross(nvec,dz4);dt4=-cross(nvec,dt4);dt3=cross(nvec,dz3);dt3=-cross(nvec,dt3);
               SurfaceIntegral(4,4)=Ae*(dt3'*dt3 +dt2'*dt2)/6 - Ae*(dt3'*dt2+ dt2'*dt3)/12;
               SurfaceIntegral(4,5)=Ae*(dt3'*dt4)/6 +Ae*(dt2'*dt2 -dt3'*dt2 -dt2'*dt4)/12;
               SurfaceIntegral(4,6)=-Ae*(dt2'*dt4)/6 +Ae*(dt2'*dt3+dt3'*dt4-dt3'*dt3)/12;
               SurfaceIntegral(5,5)=Ae*(dt4'*dt4 + dt2'*dt2)/6 -Ae*(dt4'*dt2 +dt2'*dt4)/12;
               SurfaceIntegral(5,4)=Ae*(dt4'*dt3)/6 +Ae*(dt2'*dt2 -dt2'*dt3 -dt4'*dt2)/12;
               SurfaceIntegral(5,6)=Ae*(dt2'*dt3)/6 +Ae*(dt4'*dt4 - dt4'*dt3 -dt2'*dt4)/12;
               SurfaceIntegral(6,6)=Ae*(dt3'*dt3 + dt4'*dt4)/6 -Ae*(dt4'*dt3 +dt3'*dt4)/12;
               SurfaceIntegral(6,4)=-Ae*(dt4'*dt2)/6 +Ae*(dt3'*dt2+dt4'*dt3-dt3'*dt3)/12;
               SurfaceIntegral(6,5)=Ae*(dt3'*dt2)/6 +Ae*(dt4'*dt4 - dt3'*dt4 -dt4'*dt2)/12;
        case 3,dz3=[b(3);c(3);d(3)];dz4=[b(4);c(4);d(4)];dz1=[b(1);c(1);d(1)];
               dt3=cross(nvec,dz3);dt3=-cross(nvec,dt3);dt4=cross(nvec,dz4);dt4=-cross(nvec,dt4);dt1=cross(nvec,dz1);dt1=-cross(nvec,dt1);
               SurfaceIntegral(2,2)=Ae*(dt3'*dt3 + dt1'*dt1)/6 -Ae*(dt3'*dt1 + dt1'*dt3)/12;
               SurfaceIntegral(2,3)=Ae*(dt3'*dt4)/6 + Ae*(dt1'*dt1 -dt3'*dt1 -dt1'*dt4)/12;
               SurfaceIntegral(2,6)=-Ae*(dt1'*dt4)/6 +Ae*(dt3'*dt4 +dt1'*dt3 -dt3'*dt3)/12;
               SurfaceIntegral(3,3)=Ae*(dt4'*dt4 + dt1'*dt1)/6 -Ae*(dt4'*dt1 + dt1'*dt4)/12;
               SurfaceIntegral(3,2)=Ae*(dt4'*dt3)/6 + Ae*(dt1'*dt1 -dt1'*dt3 -dt4'*dt1)/12;
               SurfaceIntegral(3,6)=Ae*(dt1'*dt3)/6+Ae*(dt4'*dt4 -dt1'*dt4 -dt4'*dt3)/12;
               SurfaceIntegral(6,6)=Ae*(dt3'*dt3 + dt4'*dt4)/6 -Ae*(dt4'*dt3 +dt3'*dt4)/12;
               SurfaceIntegral(6,2)=-Ae*(dt4'*dt1)/6 +Ae*(dt4'*dt3 +dt3'*dt1 -dt3'*dt3)/12;
               SurfaceIntegral(6,3)=Ae*(dt3'*dt1)/6+Ae*(dt4'*dt4 -dt4'*dt1 -dt3'*dt4)/12;
        case 4,dz4=[b(4);c(4);d(4)];dz2=[b(2);c(2);d(2)];dz1=[b(1);c(1);d(1)];
               dt4=cross(nvec,dz4);dt4=-cross(nvec,dt4);dt2=cross(nvec,dz2);dt2=-cross(nvec,dt2);dt1=cross(nvec,dz1);dt1=-cross(nvec,dt1);
               SurfaceIntegral(1,1)=Ae*(dt2'*dt2 + dt1'*dt1)/6 - Ae*(dt2'*dt1 + dt1'*dt2)/12;
               SurfaceIntegral(1,3)=Ae*(dt2'*dt4)/6 +Ae*(dt1'*dt1 -dt2'*dt1 -dt1'*dt4)/12;
               SurfaceIntegral(1,5)=-Ae*(dt1'*dt4)/6 +Ae*(dt2'*dt4 -dt2'*dt2 +dt1'*dt2)/12;
               SurfaceIntegral(3,3)=Ae*(dt4'*dt4 + dt1'*dt1)/6 -Ae*(dt4'*dt1 + dt1'*dt4)/12;
               SurfaceIntegral(3,1)=Ae*(dt4'*dt2)/6 +Ae*(dt1'*dt1 -dt1'*dt2 -dt4'*dt1)/12;
               SurfaceIntegral(3,5)=Ae*(dt1'*dt2)/6 +Ae*(dt4'*dt4 -dt4'*dt2 -dt1'*dt4)/12;
               SurfaceIntegral(5,5)=Ae*(dt4'*dt4 + dt2'*dt2)/6 -Ae*(dt4'*dt2 +dt2'*dt4)/12;
               SurfaceIntegral(5,1)=-Ae*(dt4'*dt1)/6 +Ae*(dt4'*dt2 -dt2'*dt2 +dt2'*dt1)/12;
               SurfaceIntegral(5,3)=Ae*(dt2'*dt1)/6 +Ae*(dt4'*dt4 -dt2'*dt4 -dt4'*dt1)/12;
    end
end
function [n] = findNormaVector(TModel,element,facetIndex),v1=TModel.Vertices(element.Vertices(1));v2=TModel.Vertices(element.Vertices(2));v3=TModel.Vertices(element.Vertices(3));v4=TModel.Vertices(element.Vertices(4));
    switch facetIndex
        case 1,vec1=[v2.X-v1.X;v2.Y-v1.Y;v2.Z-v1.Z];vec2=[v3.X-v1.X;v3.Y-v1.Y;v3.Z-v1.Z];vec3=[v1.X-v4.X;v1.Y-v4.Y;v1.Z-v4.Z];
        case 2,vec1=[v3.X-v2.X;v3.Y-v2.Y;v3.Z-v2.Z];vec2=[v4.X-v2.X;v4.Y-v2.Y;v4.Z-v2.Z];vec3=[v2.X-v1.X;v2.Y-v1.Y;v2.Z-v1.Z];
        case 3,vec1=[v3.X-v1.X;v3.Y-v1.Y;v3.Z-v1.Z];vec2=[v4.X-v1.X;v4.Y-v1.Y;v4.Z-v1.Z];vec3=[v1.X-v2.X;v1.Y-v2.Y;v1.Z-v2.Z];
        case 4,vec1=[v2.X-v1.X;v2.Y-v1.Y;v2.Z-v1.Z];vec2=[v4.X-v1.X;v4.Y-v1.Y;v4.Z-v1.Z];vec3=[v1.X-v3.X;v1.Y-v3.Y;v1.Z-v3.Z];
    end,n=cross(vec1,vec2);n=n/norm(n);if(n'*vec3>0),n=-n;end
end
function [SurfaceIntegral] = CalculateTensorIntegral(TModel,Tensor,element,facet,b,c,d,facetIndex),nvec=findNormaVector(TModel,element,facetIndex);SurfaceIntegral=zeros(6,6);Ae=facet.Surface;
    switch facetIndex
        case 1,dz1=[b(1);c(1);d(1)];dz2=[b(2);c(2);d(2)];dz3=[b(3);c(3);d(3)];
               dt1=cross(nvec,dz1);dZt1=Tensor*dt1;
               dt2=cross(nvec,dz2);dZt2=Tensor*dt2;
               dt3=cross(nvec,dz3);dZt3=Tensor*dt3;
               SurfaceIntegral(1,1)=Ae*(dt2'*dZt2 + dt1'*dZt1)/6 - Ae*(dt2'*dZt1 + dt1'*dZt2)/12;
               SurfaceIntegral(1,2)=Ae*(dt2'*dZt3)/6+Ae*(dt1'*dZt1 - dt2'*dZt1 - dt1'*dZt2)/12;
               SurfaceIntegral(1,4)=-Ae*(dt1'*dZt3)/6 +Ae*(dt1'*dZt1 - dt2'*dZt1 - dt1'*dZt3)/12;
               SurfaceIntegral(2,2)=Ae*(dt3'*dZt3 + dt1'*dZt1)/6 -Ae*(dt3'*dZt1 + dt1'*dZt3)/12;
               SurfaceIntegral(2,1)=Ae*(dt3'*dZt2)/6 + Ae*(dt1'*dZt1 -dt3'*dZt1 - dt1'*dZt2)/12;
               SurfaceIntegral(2,4)=Ae*(dt1'*dZt2)/6 + Ae*(dt3'*dZt3 - dt3'*dZt2 -dt1'*dZt3)/12;
               SurfaceIntegral(4,4)=Ae*(dt3'*dZt3 +dt2'*dZt2)/6 - Ae*(dt3'*dZt2+ dt2'*dZt3)/12;
               SurfaceIntegral(4,1)=-Ae*(dt3'*dZt1)/6 + Ae*(dt3'*dZt2 +dt2'*dZt1 -dt2'*dZt2)/12;
               SurfaceIntegral(4,2)=Ae*(dt2'*dZt1)/6 +Ae*(dt3'*dZt3 -dt3'*dZt1 -dt2'*dZt3)/12;
        case 2,dz2=[b(2);c(2);d(2)];dz4=[b(4);c(4);d(4)];dz3=[b(3);c(3);d(3)];
               dt2=cross(nvec,dz2);dZt2=Tensor*dt2;
               dt3=cross(nvec,dz3);dZt3=Tensor*dt3;
               dt4=cross(nvec,dz4);dZt4=Tensor*dt4;
               SurfaceIntegral(4,4)=Ae*(dt3'*dZt3 +dt2'*dZt2)/6 - Ae*(dt3'*dZt2+ dt2'*dZt3)/12;
               SurfaceIntegral(4,5)=Ae*(dt3'*dZt4)/6 +Ae*(dt2'*dZt2 -dt3'*dZt2 -dt2'*dZt4)/12;
               SurfaceIntegral(4,6)=-Ae*(dt2'*dZt4)/6 +Ae*(dt2'*dZt3+dt3'*dZt4-dt3'*dZt3)/12;
               SurfaceIntegral(5,5)=Ae*(dt4'*dZt4 + dt2'*dZt2)/6 -Ae*(dt4'*dZt2 +dt2'*dZt4)/12;
               SurfaceIntegral(5,4)=Ae*(dt4'*dZt3)/6 +Ae*(dt2'*dZt2 -dt2'*dZt3 -dt4'*dZt2)/12;
               SurfaceIntegral(5,6)=Ae*(dt2'*dZt3)/6 +Ae*(dt4'*dZt4 - dt4'*dZt3 -dt2'*dZt4)/12;
               SurfaceIntegral(6,6)=Ae*(dt3'*dZt3 + dt4'*dZt4)/6 -Ae*(dt4'*dZt3 +dt3'*dZt4)/12;
               SurfaceIntegral(6,4)=-Ae*(dt4'*dZt2)/6 +Ae*(dt3'*dZt2+dt4'*dZt3-dt3'*dZt3)/12;
               SurfaceIntegral(6,5)=Ae*(dt3'*dZt2)/6 +Ae*(dt4'*dZt4 - dt3'*dZt4 -dt4'*dZt2)/12;
        case 3,dz3=[b(3);c(3);d(3)];dz4=[b(4);c(4);d(4)];dz1=[b(1);c(1);d(1)];
               dt1=cross(nvec,dz1);dZt1=Tensor*dt1;
               dt4=cross(nvec,dz4);dZt4=Tensor*dt4;
               dt3=cross(nvec,dz3);dZt3=Tensor*dt3;
               SurfaceIntegral(2,2)=Ae*(dt3'*dZt3 + dt1'*dZt1)/6 -Ae*(dt3'*dZt1 + dt1'*dZt3)/12;
               SurfaceIntegral(2,3)=Ae*(dt3'*dZt4)/6 + Ae*(dt1'*dZt1 -dt3'*dZt1 -dt1'*dZt4)/12;
               SurfaceIntegral(2,6)=-Ae*(dt1'*dZt4)/6 +Ae*(dt3'*dZt4 +dt1'*dZt3 -dt3'*dZt3)/12;
               SurfaceIntegral(3,3)=Ae*(dt4'*dZt4 + dt1'*dZt1)/6 -Ae*(dt4'*dZt1 + dt1'*dZt4)/12;
               SurfaceIntegral(3,2)=Ae*(dt4'*dZt3)/6 + Ae*(dt1'*dZt1 -dt1'*dZt3 -dt4'*dZt1)/12;
               SurfaceIntegral(3,6)=Ae*(dt1'*dZt3)/6+Ae*(dt4'*dZt4 -dt1'*dZt4 -dt4'*dZt3)/12;
               SurfaceIntegral(6,6)=Ae*(dt3'*dZt3 + dt4'*dZt4)/6 -Ae*(dt4'*dZt3 +dt3'*dZt4)/12;
               SurfaceIntegral(6,2)=-Ae*(dt4'*dZt1)/6 +Ae*(dt4'*dZt3 +dt3'*dZt1 -dt3'*dZt3)/12;
               SurfaceIntegral(6,3)=Ae*(dt3'*dZt1)/6+Ae*(dt4'*dZt4 -dt4'*dZt1 -dt3'*dZt4)/12;
        case 4,dz4=[b(4);c(4);d(4)];dz2=[b(2);c(2);d(2)];dz1=[b(1);c(1);d(1)];
               dt4=cross(nvec,dz4);dZt4=Tensor*dt4;
               dt2=cross(nvec,dz2);dZt2=Tensor*dt2;
               dt1=cross(nvec,dz1);dZt1=Tensor*dt1;
               SurfaceIntegral(1,1)=Ae*(dt2'*dZt2 + dt1'*dZt1)/6 - Ae*(dt2'*dZt1 + dt1'*dZt2)/12;
               SurfaceIntegral(1,3)=Ae*(dt2'*dZt4)/6 +Ae*(dt1'*dZt1 -dt2'*dZt1 -dt1'*dZt4)/12;
               SurfaceIntegral(1,5)=-Ae*(dt1'*dZt4)/6 +Ae*(dt2'*dZt4 -dt2'*dZt2 +dt1'*dZt2)/12;
               SurfaceIntegral(3,3)=Ae*(dt4'*dZt4 + dt1'*dZt1)/6 -Ae*(dt4'*dZt1 + dt1'*dZt4)/12;
               SurfaceIntegral(3,1)=Ae*(dt4'*dZt2)/6 +Ae*(dt1'*dZt1 -dt1'*dZt2 -dt4'*dZt1)/12;
               SurfaceIntegral(3,5)=Ae*(dt1'*dZt2)/6 +Ae*(dt4'*dZt4 -dt4'*dZt2 -dt1'*dZt4)/12;
               SurfaceIntegral(5,5)=Ae*(dt4'*dZt4 + dt2'*dZt2)/6 -Ae*(dt4'*dZt2 +dt2'*dZt4)/12;
               SurfaceIntegral(5,1)=-Ae*(dt4'*dZt1)/6 +Ae*(dt4'*dZt2 -dt2'*dZt2 +dt2'*dZt1)/12;
               SurfaceIntegral(5,3)=Ae*(dt2'*dZt1)/6 +Ae*(dt4'*dZt4 -dt2'*dZt4 -dt4'*dZt1)/12;
    end,SurfaceIntegral=-SurfaceIntegral;
end
