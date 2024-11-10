%==========================================================================
%{
          Finite Element Assembly Function for E-B (Maxwell Equations) 
                        Excitation Formulation

   i.Isotropic and Anisotropic Sub Functions :
   Zhu, Y. and Cangellaris, A.C. eds., 2006. Multigrid finite element methods 
   for electromagnetic field modeling (Vol. 28). John Wiley & Sons.
   
   Ανάπτυξη Υπολογιστικών Τεχνικών με τη Χρήση Διατυπώσεων Μεικτών Πεπερασμένων 
   Στοιχείων στο Σύστημα Εξισώσεων του Maxwell για Προβλήματα Ηλεκτρομαγνητικής 
   Διάδοσης, Σαλονικιός Βασίλειος,2024 (Ph.D. Thesis)

   ii. Bianisotropic Media 

   Nitas, M., Salonikios, V., Amanatiadis, S. and Arslanagic, S., 2023, March. 
   Field-Flux Finite Element Formulation for Wave Propagation in Bianisotropic 
   Media. In 2023 17th European Conference on Antennas and Propagation (EuCAP) 
   (pp. 1-4). IEEE.
    
%}
%==========================================================================
function [TModel] = Excitation_Assembly(TModel),Assembly=TModel.Assembled;Frequency=TModel.Frequency;
%--------------------------- Single Frequency -----------------------------
    if(Frequency.NF==1),frequency=Frequency.Frequency;
    %------------------------- Dirichlet Boundary Excitation --------------
    if(Assembly.IsDir)
        %--------------------- Sparse Matrix Initializations --------------
        if(Assembly.Bian),N=400;
            IA=zeros(N*TModel.NumberOfElements,1);JA=IA;TE=IA;TB=IA;F=IA;A=IA;counterA=0;
            P=IA;TA=IA;TC=IA;TS=IA;TP=IA;TG=IA;TBC=IA;
            IB=zeros(N*TModel.NumberOfElements,1);JB=IB;TEB=IB;TBB=IB;FB=IB;AB=IB;counterB=0;
            PB=IB;TAB=IB;TCB=IB;
        else,N=400;
            IA=zeros(N*TModel.NumberOfElements,1);JA=IA;TE=IA;TB=IA;F=IA;A=IA;counterA=0;
            IB=zeros(N*TModel.NumberOfElements,1);JB=IB;TEB=IB;TBB=IB;FB=IB;AB=IB;counterB=0;
            TS=IA;TP=IA;TG=IA;TBC=IA;
        end
        %------------------------- Domains Assembly -----------------------
        for id = 1 :numel (TModel.Domains),domain = TModel.Domains(id);medium=domain.Medium;
            switch medium.Type
                case "Iso" ,[IA,JA,TE,TB,F,A,TS,TP,TG,TBC,IB,JB,TEB,TBB,FB,AB,counterA,counterB] = DirichletExcitation_Iso(IA,JA,TE,TB,F,A,TS,TP,TG,TBC,IB,JB,TEB,TBB,FB,AB,counterA,counterB,TModel,domain);
                case "Anis",[IA,JA,TE,TB,F,A,TS,TP,TG,TBC,IB,JB,TEB,TBB,FB,AB,counterA,counterB] = DirichletExcitation_Ani(IA,JA,TE,TB,F,A,TS,TP,TG,TBC,IB,JB,TEB,TBB,FB,AB,counterA,counterB,TModel,domain);
                case "Bian",[IA,JA,TE,TB,F,A,TS,TP,TG,TBC,P,TA,TC,IB,JB,TEB,TBB,FB,AB,PB,TAB,TCB,counterA,counterB] = DirichletExcitation_Bia(IA,JA,TE,TB,F,A,TS,TP,TG,TBC,P,TA,TC,IB,JB,TEB,TBB,FB,AB,PB,TAB,TCB,counterA,counterB,TModel,domain);
            end
        end
        IB(counterB)=max([TModel.Facets.UknownIndex]);JB(counterB)=max([TModel.Facets.KnownIndex]);
        %------------------------- Sparse Matrix Construction -------------
        nonZeros=nnz(IA); IA=IA(1:nonZeros);JA=JA(1:nonZeros);        nonZeros_B=nnz(IB);IB=IB(1:nonZeros_B);JB=JB(1:nonZeros_B);
       
        TE=TE(1:nonZeros);      TS=TS(1:nonZeros);          TEB=TEB(1:nonZeros_B);
        TB=TB(1:nonZeros);      TP=TP(1:nonZeros);          TBB=TBB(1:nonZeros_B);
        F=F(1:nonZeros);        TG=TG(1:nonZeros);          FB=FB(1:nonZeros_B);
        A=A(1:nonZeros);        TBC=TBC(1:nonZeros);        AB=AB(1:nonZeros_B);
      
        TE=sparse(IA,JA,TE);    TS=sparse(IA,JA,TS);        TEB=sparse(IB,JB,TEB);
        TB=sparse(IA,JA,TB);    TP=sparse(IA,JA,TP);        TBB=sparse(IB,JB,TBB);
        A=sparse(IA,JA,A);      TG=sparse(IA,JA,TG);        AB=sparse(IB,JB,AB);
        F=sparse(IA,JA,F);      TBC=sparse(IA,JA,TBC);      FB=sparse(IB,JB,FB);
        
        if(~Assembly.Bian),Assembly=Assembly.DirichletExcitation(TE,TB,A,F,TS,TP,TG,TBC,TEB,TBB,AB,FB,TModel,frequency);TModel.Assembled=Assembly;
        else,P=P(1:nonZeros);       PB=PB(1:nonZeros_B);
             TA=TA(1:nonZeros);     TAB=TAB(1:nonZeros_B);
             TC=TC(1:nonZeros);     TCB=TCB(1:nonZeros_B);

             P=sparse(IA,JA,P);     PB=sparse(IB,JB,PB);
             TA=sparse(IA,JA,TA);   TAB=sparse(IB,JB,TAB);
             TC=sparse(IA,JA,TC);   TCB=sparse(IB,JB,TCB);

             Assembly=Assembly.BianisotropicDirichletExcitation(TE,TB,A,F,TS,TP,TG,TBC,P,TA,TC,TEB,TBB,AB,FB,PB,TAB,TCB,TModel,frequency);TModel.Assembled=Assembly;
        end
    else
     %----------------------------- Port Boundary Excitation -------------- 
        if(Assembly.Bian),N=400;
            IA=zeros(N*TModel.NumberOfElements,1);JA=IA;TE=IA;TB=IA;F=IA;A=IA;counterA=0;
            P=IA;TA=IA;TC=IA;TS=IA;TP=IA;TG=IA;TBC=IA;
        else,N=400;
            IA=zeros(N*TModel.NumberOfElements,1);JA=IA;TE=IA;TB=IA;F=IA;A=IA;counterA=0;
            TS=IA;TP=IA;TG=IA;TBC=IA;
        end
    end
        %------------------------- Domains Assembly -----------------------
         for id = 1 :numel (TModel.Domains),domain = TModel.Domains(id);medium=domain.Medium;
            switch medium.Type
                case "Iso" ,[IA,JA,TE,TB,F,A,TS,TP,TG,TBC,counterA] = PortExcitation_Iso(IA,JA,TE,TB,F,A,TS,TP,TG,TBC,counterA,TModel,domain);
                case "Anis",[IA,JA,TE,TB,F,A,TS,TP,TG,TBC,counterA] = PortExcitation_Ani(IA,JA,TE,TB,F,A,TS,TP,TG,TBC,counterA,TModel,domain);
                case "Bian",[IA,JA,TE,TB,F,A,TS,TP,TG,TBC,P,TA,TC,counterA] = PortExcitation_Bia(IA,JA,TE,TB,F,A,TS,TP,TG,TBC,P,TA,TC,counterA,TModel,domain);
            end
         end
         %-----------------------------------------------------------------
        nonZeros=nnz(IA);IA=IA(1:nonZeros);JA=JA(1:nonZeros);
        TE=TE(1:nonZeros);      TS=TS(1:nonZeros);
        TB=TB(1:nonZeros);      TP=TP(1:nonZeros);
        F=F(1:nonZeros);        TG=TG(1:nonZeros);
        A=A(1:nonZeros);        TBC=TBC(1:nonZeros);

        TE=sparse(IA,JA,TE);        TS=sparse(IA,JA,TS);
        TB=sparse(IA,JA,TB);        TP=sparse(IA,JA,TP);
        A=sparse(IA,JA,A);          TG=sparse(IA,JA,TG);
        F=sparse(IA,JA,F);          TBC=sparse(IA,JA,TBC);
        if(~Assembly.Bian),Assembly=Assembly.PortExcitation(TE,TB,A,F,TS,TP,TG,TBC,TModel,frequency);TModel.Assembled=Assembly;
        else
             P=P(1:nonZeros);       P=sparse(IA,JA,P);
             TA=TA(1:nonZeros);     TA=sparse(IA,JA,TA);
             TC=TC(1:nonZeros);     TC=sparse(IA,JA,TC);
             Assembly=Assembly.BianisotropicPortExcitation(TE,TB,A,F,TS,TP,TG,TBC,P,TA,TC,TModel,frequency);TModel.Assembled=Assembly;
        end
%---------------------------- Multi Frequency -----------------------------
    else
        for ii=1:Frequency.NF,frequency=Frequency.Frequency(ii);
            if(Assembly.IsDir)
            %------------------------- Dirichlet Boundary Excitation ------
                if(Assembly.Bian),N=400;
                    IA=zeros(N*TModel.NumberOfElements,1);JA=IA;TE=IA;TB=IA;F=IA;A=IA;counterA=0;
                    P=IA;TA=IA;TC=IA;TS=IA;TP=IA;TG=IA;TBC=IA;
                    IB=zeros(N*TModel.NumberOfElements,1);JB=IB;TEB=IB;TBB=IB;FB=IB;AB=IB;counterB=0;
                    PB=IB;TAB=IB;TCB=IB;
                else,N=400;
                    IA=zeros(N*TModel.NumberOfElements,1);JA=IA;TE=IA;TB=IA;F=IA;A=IA;counterA=0;
                    IB=zeros(N*TModel.NumberOfElements,1);JB=IB;TEB=IB;TBB=IB;FB=IB;AB=IB;counterB=0;
                    TS=IA;TP=IA;TG=IA;TBC=IA;
                end
                 %------------------------- Domains Assembly --------------
                for id = 1 :numel (TModel.Domains),domain = TModel.Domains(id);medium=domain.Medium;
                    switch medium.Type
                        case "Iso" ,[IA,JA,TE,TB,F,A,TS,TP,TG,TBC,IB,JB,TEB,TBB,FB,AB,counterA,counterB] = DirichletExcitation_Iso(IA,JA,TE,TB,F,A,TS,TP,TG,TBC,IB,JB,TEB,TBB,FB,AB,counterA,counterB,TModel,domain,ii);
                        case "Anis",[IA,JA,TE,TB,F,A,TS,TP,TG,TBC,IB,JB,TEB,TBB,FB,AB,counterA,counterB] = DirichletExcitation_Ani(IA,JA,TE,TB,F,A,TS,TP,TG,TBC,IB,JB,TEB,TBB,FB,AB,counterA,counterB,TModel,domain,ii);
                        case "Bian",[IA,JA,TE,TB,F,A,TS,TP,TG,TBC,P,TA,TC,IB,JB,TEB,TBB,FB,AB,PB,TAB,TCB,counterA,counterB] = DirichletExcitation_Bia(IA,JA,TE,TB,F,A,TS,TP,TG,TBC,P,TA,TC,IB,JB,TEB,TBB,FB,AB,PB,TAB,TCB,counterA,counterB,TModel,domain,ii);
                    end
                end
                 %------------------------- Sparse Matrix Construction -------------
                nonZeros=nnz(IA);IA=IA(1:nonZeros);JA=JA(1:nonZeros);nonZeros_B=nnz(IB);IB=IB(1:nonZeros_B);JB=JB(1:nonZeros_B);
                TE=TE(1:nonZeros);      TS=TS(1:nonZeros);      TEB=TEB(1:nonZeros_B);
                TB=TB(1:nonZeros);      TP=TP(1:nonZeros);      TBB=TBB(1:nonZeros_B);
                F=F(1:nonZeros);        TG=TG(1:nonZeros);      FB=FB(1:nonZeros_B);
                A=A(1:nonZeros);        TBC=TBC(1:nonZeros);    AB=AB(1:nonZeros_B);
                
                TE=sparse(IA,JA,TE);        TS=sparse(IA,JA,TS);        TEB=sparse(IB,JB,TEB);
                TB=sparse(IA,JA,TB);        TP=sparse(IA,JA,TP);        TBB=sparse(IB,JB,TBB);
                A=sparse(IA,JA,A);          TG=sparse(IA,JA,TG);        AB=sparse(IB,JB,AB);
                F=sparse(IA,JA,F);          TBC=sparse(IA,JA,TBC);      FB=sparse(IB,JB,FB);
               
                if(~Assembly.Bian),Assembly=Assembly.DirichletExcitation(TE,TB,A,F,TS,TP,TG,TBC,TEB,TBB,AB,FB,TModel,frequency,ii);TModel.Assembled=Assembly;
                        else,P=P(1:nonZeros);       PB=PB(1:nonZeros_B);
                             TA=TA(1:nonZeros);     TAB=TAB(1:nonZeros_B);
                             TC=TC(1:nonZeros);     TCB=TCB(1:nonZeros_B);
                             
                             P=sparse(IA,JA,P);     PB=sparse(IB,JB,PB);
                             TA=sparse(IA,JA,TA);   TAB=sparse(IB,JB,TAB);
                             TC=sparse(IA,JA,TC);   TCB=sparse(IB,JB,TCB);
                             Assembly=Assembly.BianisotropicDirichletExcitation(TE,TB,A,F,TS,TP,TG,TBC,P,TA,TC,TEB,TBB,AB,FB,PB,TAB,TCB,TModel,frequency,ii);TModel.Assembled=Assembly;
                end
            else
                %----------------------------- Port Boundary Excitation -------------- 
                if(Assembly.Bian),N=400;
                    IA=zeros(N*TModel.NumberOfElements,1);JA=IA;TE=IA;TB=IA;F=IA;A=IA;counterA=0;
                    P=IA;TA=IA;TC=IA;TS=IA;TP=IA;TG=IA;TBC=IA;
                else,N=400;
                    IA=zeros(N*TModel.NumberOfElements,1);JA=IA;TE=IA;TB=IA;F=IA;A=IA;counterA=0;
                    TS=IA;TP=IA;TG=IA;TBC=IA;
                end
                %------------------------- Domains Assembly ---------------
                for id = 1 :numel (TModel.Domains),domain = TModel.Domains(id);medium=domain.Medium;
                      switch medium.Type
                            case "Iso" ,[IA,JA,TE,TB,F,A,TS,TP,TG,TBC,counterA] = PortExcitation_Iso(IA,JA,TE,TB,F,A,TS,TP,TG,TBC,counterA,TModel,domain,ii);
                            case "Anis",[IA,JA,TE,TB,F,A,TS,TP,TG,TBC,counterA] = PortExcitation_Ani(IA,JA,TE,TB,F,A,TS,TP,TG,TBC,counterA,TModel,domain,ii);
                            case "Bian",[IA,JA,TE,TB,F,A,TS,TP,TG,TBC,P,TA,TC,counterA] = PortExcitation_Bia(IA,JA,TE,TB,F,A,TS,TP,TG,TBC,P,TA,TC,counterA,TModel,domain,ii);
                      end
                end
                %-----------------------------------------------------------
                nonZeros=nnz(IA);IA=IA(1:nonZeros);JA=JA(1:nonZeros);
                TE=TE(1:nonZeros);      TS=TS(1:nonZeros);
                TB=TB(1:nonZeros);      TP=TP(1:nonZeros);
                F=F(1:nonZeros);        TG=TG(1:nonZeros);
                A=A(1:nonZeros);        TBC=TBC(1:nonZeros);

                TE=sparse(IA,JA,TE);        TS=sparse(IA,JA,TS);
                TB=sparse(IA,JA,TB);        TP=sparse(IA,JA,TP);
                A=sparse(IA,JA,A);          TG=sparse(IA,JA,TG);
                F=sparse(IA,JA,F);          TBC=sparse(IA,JA,TBC);
                if(~Assembly.Bian),Assembly=Assembly.PortExcitation(TE,TB,A,F,TS,TP,TG,TBC,TModel,frequency,ii);TModel.Assembled=Assembly;
                else,P=P(1:nonZeros);       P=sparse(IA,JA,P);
                     TA=TA(1:nonZeros);     TA=sparse(IA,JA,TA);
                     TC=TC(1:nonZeros);     TC=sparse(IA,JA,TC);
                     Assembly=Assembly.BianisotropicPortExcitation(TE,TB,A,F,TS,TP,TG,TBC,P,TA,TC,TModel,frequency,ii);TModel.Assembled=Assembly;
                end
           end
        end
    end
end
%====================== Dirichlet Excitation Assembly =====================
function [IA,JA,TE,TB,FF,AA,TS,TP,TG,TBC,IB,JB,TEB,TBB,FB,AB,counterA,counterB] =  DirichletExcitation_Iso(varargin),ElectromagneticConstants;GaussianQuadratture3D;
    if(nargin==20)
        IA=varargin{1};     TS=varargin{7};         IB=varargin{11};        counterA=varargin{17};
        JA=varargin{2};     TP=varargin{8};         JB=varargin{12};        counterB=varargin{18};
        TE=varargin{3};     TG=varargin{9};         TEB=varargin{13};       TModel=varargin{19};
        TB=varargin{4};     TBC=varargin{10};       TBB=varargin{14};       domain=varargin{20};
        FF=varargin{5};                             FB=varargin{15};
        AA=varargin{6};                             AB=varargin{16};
        medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;freq=TModel.Frequency.Frequency;
    elseif(nargin==21)
        IA=varargin{1};     TS=varargin{7};         IB=varargin{11};        counterA=varargin{17};
        JA=varargin{2};     TP=varargin{8};         JB=varargin{12};        counterB=varargin{18};
        TE=varargin{3};     TG=varargin{9};         TEB=varargin{13};       TModel=varargin{19};
        TB=varargin{4};     TBC=varargin{10};       TBB=varargin{14};       domain=varargin{20};
        FF=varargin{5};                             FB=varargin{15};        freqIndex=varargin{21};
        AA=varargin{6};                             AB=varargin{16};
        medium=domain.Medium;freq=TModel.Frequency.Frequency(freqIndex);
        if(medium.IsDispersive),epsilon=medium.Epsilon(freqIndex);mu=medium.Mu(freqIndex);imu=mu^-1;
        else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;
        end
    end,k0=2*pi*freq/c0;
    for ie=1:numel(domain.Elements),element=TModel.Elements(domain.Elements(ie));
        b=element.Bs;c=element.Cs;d=element.Ds;Ve=element.Volume;edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];
        edgeSigns=[element.EdgeSigns];facetSigns=[element.FacetSigns];edgeLengths=[edges.Length];facetSurfaces=[facets.Surface];
        Te=zeros(6,6);Tb=zeros(4,4);Am=zeros(6,4);Fm=zeros(4,6);Tp=zeros(6,6);Ts=zeros(6,6);Tbc=zeros(6,6);Tg=zeros(6,6);
        for kt =1: numel(Weights),zeta(1)=Points(1,kt);zeta(2)=Points(2,kt);zeta(3)=Points(3,kt);zeta(4)=Points(4,kt);
                 %---------------------- Edge Basis Functions -------------
                 wx(1)=zeta(1)*b(2)-zeta(2)*b(1);       wy(1)=zeta(1)*c(2)-zeta(2)*c(1);         wz(1)=zeta(1)*d(2)-zeta(2)*d(1);  
                 wx(2)=zeta(1)*b(3)-zeta(3)*b(1);       wy(2)=zeta(1)*c(3)-zeta(3)*c(1);         wz(2)=zeta(1)*d(3)-zeta(3)*d(1);
                 wx(3)=zeta(1)*b(4)-zeta(4)*b(1);       wy(3)=zeta(1)*c(4)-zeta(4)*c(1);         wz(3)=zeta(1)*d(4)-zeta(4)*d(1);
                 wx(4)=zeta(2)*b(3)-zeta(3)*b(2);       wy(4)=zeta(2)*c(3)-zeta(3)*c(2);         wz(4)=zeta(2)*d(3)-zeta(3)*d(2);
                 wx(5)=zeta(2)*b(4)-zeta(4)*b(2);       wy(5)=zeta(2)*c(4)-zeta(4)*c(2);         wz(5)=zeta(2)*d(4)-zeta(4)*d(2);
                 wx(6)=zeta(3)*b(4)-zeta(4)*b(3);       wy(6)=zeta(3)*c(4)-zeta(4)*c(3);         wz(6)=zeta(3)*d(4)-zeta(4)*d(3);
                 %------------------ Edge Basis Functions Rotations -------
                 rwx(1)=2*(c(1)*d(2)-d(1)*c(2));        rwy(1)=2*(d(1)*b(2)-b(1)*d(2));          rwz(1)=2*(b(1)*c(2)-c(1)*b(2));
                 rwx(2)=2*(c(1)*d(3)-d(1)*c(3));        rwy(2)=2*(d(1)*b(3)-b(1)*d(3));          rwz(2)=2*(b(1)*c(3)-c(1)*b(3));
                 rwx(3)=2*(c(1)*d(4)-d(1)*c(4));        rwy(3)=2*(d(1)*b(4)-b(1)*d(4));          rwz(3)=2*(b(1)*c(4)-c(1)*b(4));
                 rwx(4)=2*(c(2)*d(3)-d(2)*c(3));        rwy(4)=2*(d(2)*b(3)-b(2)*d(3));          rwz(4)=2*(b(2)*c(3)-c(2)*b(3));
                 rwx(5)=2*(c(2)*d(4)-d(2)*c(4));        rwy(5)=2*(d(2)*b(4)-b(2)*d(4));          rwz(5)=2*(b(2)*c(4)-c(2)*b(4));
                 rwx(6)=2*(c(3)*d(4)-d(3)*c(4));        rwy(6)=2*(d(3)*b(4)-b(3)*d(4));          rwz(6)=2*(b(3)*c(4)-c(3)*b(4));
                 %---------------------- Facet Basis Functions ------------
                 wfx(1)=2*zeta(3)*(c(1)*d(2)-c(2)*d(1))+2*zeta(1)*(c(2)*d(3)-c(3)*d(2))+2*zeta(2)*(c(3)*d(1)-c(1)*d(3));    wfy(1)=2*zeta(3)*(d(1)*b(2)-d(2)*b(1))+2*zeta(1)*(d(2)*b(3)-d(3)*b(2))+2*zeta(2)*(d(3)*b(1)-d(1)*b(3));     wfz(1)=2*zeta(3)*(b(1)*c(2)-b(2)*c(1))+2*zeta(1)*(b(2)*c(3)-b(3)*c(2))+2*zeta(2)*(b(3)*c(1)-b(1)*c(3));
                 wfx(2)=2*zeta(3)*(c(2)*d(4)-c(4)*d(2))+2*zeta(2)*(c(4)*d(3)-c(3)*d(4))+2*zeta(4)*(c(3)*d(2)-c(2)*d(3));    wfy(2)=2*zeta(3)*(d(2)*b(4)-d(4)*b(2))+2*zeta(2)*(d(4)*b(3)-d(3)*b(4))+2*zeta(4)*(d(3)*b(2)-d(2)*b(3));     wfz(2)=2*zeta(3)*(b(2)*c(4)-b(4)*c(2))+2*zeta(2)*(b(4)*c(3)-b(3)*c(4))+2*zeta(4)*(b(3)*c(2)-b(2)*c(3));
                 wfx(3)=2*zeta(1)*(c(3)*d(4)-c(4)*d(3))+2*zeta(3)*(c(4)*d(1)-c(1)*d(4))+2*zeta(4)*(c(1)*d(3)-c(3)*d(1));    wfy(3)=2*zeta(1)*(d(3)*b(4)-d(4)*b(3))+2*zeta(3)*(d(4)*b(1)-d(1)*b(4))+2*zeta(4)*(d(1)*b(3)-d(3)*b(1));     wfz(3)=2*zeta(1)*(b(3)*c(4)-b(4)*c(3))+2*zeta(3)*(b(4)*c(1)-b(1)*c(4))+2*zeta(4)*(b(1)*c(3)-b(3)*c(1));
                 wfx(4)=2*zeta(1)*(c(4)*d(2)-c(2)*d(4))+2*zeta(4)*(c(2)*d(1)-c(1)*d(2))+2*zeta(2)*(c(1)*d(4)-c(4)*d(1));    wfy(4)=2*zeta(1)*(d(4)*b(2)-d(2)*b(4))+2*zeta(4)*(d(2)*b(1)-d(1)*b(2))+2*zeta(2)*(d(1)*b(4)-d(4)*b(1));     wfz(4)=2*zeta(1)*(b(4)*c(2)-b(2)*c(4))+2*zeta(4)*(b(2)*c(1)-b(1)*c(2))+2*zeta(2)*(b(1)*c(4)-b(4)*c(1));
                 %===================== Basis Function Scaling ============
                 switch TModel.Assembled.E_Scaling
                     case 1,wx=wx.*edgeLengths;     rwx=rwx.*edgeLengths;
                            wy=wy.*edgeLengths;     rwy=rwy.*edgeLengths;
                            wz=wz.*edgeLengths;     rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch TModel.Assembled.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*TModel.Assembled.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*TModel.Assembled.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*TModel.TModel.Assembled.ScalingVec_B(facets(ii).Index);end
                 end
                 %==================== Local Matrices =====================
                 for ii=1:6,for jj=1:6,Te(ii,jj)=Te(ii,jj)+Weights(kt)*epsilon*(wx(ii)*wx(jj)+wy(ii)*wy(jj)+wz(ii)*wz(jj));end,end
                 for ii=1:4,for jj=1:4,Tb(ii,jj)=Tb(ii,jj)+Weights(kt)*(wfx(ii)*wfx(jj)+wfy(ii)*wfy(jj)+wfz(ii)*wfz(jj));end,end
                 for ii=1:4,for jj=1:6,Fm(ii,jj)=Fm(ii,jj)+Weights(kt)*(wfx(ii)*rwx(jj)+wfy(ii)*rwy(jj)+wfz(ii)*rwz(jj));end,end
                 for ii=1:6,for jj=1:4,Am(ii,jj)=Am(ii,jj)+Weights(kt)*imu*(rwx(ii)*wfx(jj)+rwy(ii)*wfy(jj)+rwz(ii)*wfz(jj));end,end
        end,Te=Te*Ve;Tb=Tb*Ve;Am=Am*Ve;Fm=Fm*Ve;
        %======================== Surface Integral Terms ==================
        for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                switch boundary.Type
                    case "ABC"
                        Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(medium.IsDispersive),ZW=medium.WaveImpedance(freqIndex);else,ZW=medium.WaveImpedance;end,Ts=(ZW^-1)*Ts;Ts=Ts*Z0;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "ABB"
                        Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(boundary.Dispersive),beta=boundary.Param(freqIndex);else,beta=boundary.Param;end,k0=2*pi*freq/c0;Ts=Ts*beta/k0;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "ABZ"
                       if(boundary.Tensor)
                            if(boundary.Dispersive),Z=boundary.Param{freqIndex};
                            else,Z=boundary.Param;
                            end
                            Ts=CalculateTensorIntegral(TModel,Z^-1,element,facets(kk),b,c,d,kk);Ts=Ts*Z0;
                        else
                            if(boundary.Dispersive),Z=boundary.Param(freqIndex);
                            else,Z=boundary.Param;
                            end
                            Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=Ts*Z^-1;Ts=Ts*Z0;
                        end
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "GRA"
                        Tg = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(boundary.Dispersive),cond=boundary.Param(freqIndex);else,cond=boundary.Param;end,Tg=Tg.*cond;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "IBC"
                        if(boundary.Dispersive),cond=boundary.Param(freqIndex);else,cond=boundary.Param;end,Tbc=Tbc.*cond;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "POR"
                        if(boundary.Dispersive)
                            if(boundary.PortParamType==0)
                                Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Ts=Ts*boundary.Param(freqIndex)/k0;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==1),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Zp=boundary.Param(freqIndex);
                                Ts=(Z0/Zp)*Ts;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==2)
                                Zp=boundary.Param{freqIndex};
                                Ts=CalculateTensorIntegral(TModel,Zp^-1,element,facets(kk),b,c,d,kk);
                                Ts=Z0*Ts;Tp=2*1i*Ts;
                            end
                        else
                            if(boundary.PortParamType==0)
                                Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Ts=boundary.Param*Ts/k0;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==1),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Zp=boundary.Param; 
                                Ts=(Z0/Zp)*Ts;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==2),Zp=boundary.Param;Ts=CalculateTensorIntegral(TModel,Zp^-1,element,facets(kk),b,c,d,kk);
                                Ts=Z0*Ts;Tp=2*1i*Ts;
                            end
                        end
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6
                                        for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);
                                                   Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);
                                        end
                                    end
                        end
                end
            end
        end
        %==================== Global Matrix Assembly ======================
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TP(counterA)=dsi*dsj*si*sj*Tp(ii,jj);
                end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;TBB(counterB)=si*sj*Te(ii,jj);end
            end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;AB(counterB)=si*sj*Am(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;FB(counterB)=si*sj*Fm(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj); end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;TBB(counterB)=si*sj*dsi*dsj*Tb(ii,jj);end
            end
        end
    end
end
%--------------------------------------------------------------------------
function [IA,JA,TE,TB,FF,AA,TS,TP,TG,TBC,IB,JB,TEB,TBB,FB,AB,counterA,counterB] =  DirichletExcitation_Ani(varargin),ElectromagneticConstants;GaussianQuadratture3D;
    if(nargin==20)
        IA=varargin{1};     TS=varargin{7};         IB=varargin{11};        counterA=varargin{17};
        JA=varargin{2};     TP=varargin{8};         JB=varargin{12};        counterB=varargin{18};
        TE=varargin{3};     TG=varargin{9};         TEB=varargin{13};       TModel=varargin{19};
        TB=varargin{4};     TBC=varargin{10};       TBB=varargin{14};       domain=varargin{20};
        FF=varargin{5};                             FB=varargin{15};
        AA=varargin{6};                             AB=varargin{16};
        medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;freq=TModel.Frequency.Frequency;
    elseif(nargin==21)
        IA=varargin{1};     TS=varargin{7};         IB=varargin{11};        counterA=varargin{17};
        JA=varargin{2};     TP=varargin{8};         JB=varargin{12};        counterB=varargin{18};
        TE=varargin{3};     TG=varargin{9};         TEB=varargin{13};       TModel=varargin{19};
        TB=varargin{4};     TBC=varargin{10};       TBB=varargin{14};       domain=varargin{20};
        FF=varargin{5};                             FB=varargin{15};        freqIndex=varargin{21};
        AA=varargin{6};                             AB=varargin{16};
        medium=domain.Medium;freq=TModel.Frequency.Frequency(freqIndex);
        if(medium.IsDispersive),epsilon=medium.Epsilon{freqIndex};mu=medium.Mu{freqIndex};imu=mu^-1;
        else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;
        end
    end,k0=2*pi*freq/c0;
    for ie=1:numel(domain.Elements),element=TModel.Elements(domain.Elements(ie));
        b=element.Bs;c=element.Cs;d=element.Ds;Ve=element.Volume;edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];
        edgeSigns=[element.EdgeSigns];facetSigns=[element.FacetSigns];edgeLengths=[edges.Length];facetSurfaces=[facets.Surface];
        Te=zeros(6,6);Tb=zeros(4,4);Am=zeros(6,4);Fm=zeros(4,6);Tp=zeros(6,6);Ts=zeros(6,6);Tbc=zeros(6,6);Tg=zeros(6,6);
        for kt =1: numel(Weights),zeta(1)=Points(1,kt);zeta(2)=Points(2,kt);zeta(3)=Points(3,kt);zeta(4)=Points(4,kt);
                 %---------------------- Edge Basis Functions -------------
                 wx(1)=zeta(1)*b(2)-zeta(2)*b(1);       wy(1)=zeta(1)*c(2)-zeta(2)*c(1);         wz(1)=zeta(1)*d(2)-zeta(2)*d(1);  
                 wx(2)=zeta(1)*b(3)-zeta(3)*b(1);       wy(2)=zeta(1)*c(3)-zeta(3)*c(1);         wz(2)=zeta(1)*d(3)-zeta(3)*d(1);
                 wx(3)=zeta(1)*b(4)-zeta(4)*b(1);       wy(3)=zeta(1)*c(4)-zeta(4)*c(1);         wz(3)=zeta(1)*d(4)-zeta(4)*d(1);
                 wx(4)=zeta(2)*b(3)-zeta(3)*b(2);       wy(4)=zeta(2)*c(3)-zeta(3)*c(2);         wz(4)=zeta(2)*d(3)-zeta(3)*d(2);
                 wx(5)=zeta(2)*b(4)-zeta(4)*b(2);       wy(5)=zeta(2)*c(4)-zeta(4)*c(2);         wz(5)=zeta(2)*d(4)-zeta(4)*d(2);
                 wx(6)=zeta(3)*b(4)-zeta(4)*b(3);       wy(6)=zeta(3)*c(4)-zeta(4)*c(3);         wz(6)=zeta(3)*d(4)-zeta(4)*d(3);
                 %------------------ Edge Basis Functions Rotations -------
                 rwx(1)=2*(c(1)*d(2)-d(1)*c(2));        rwy(1)=2*(d(1)*b(2)-b(1)*d(2));          rwz(1)=2*(b(1)*c(2)-c(1)*b(2));
                 rwx(2)=2*(c(1)*d(3)-d(1)*c(3));        rwy(2)=2*(d(1)*b(3)-b(1)*d(3));          rwz(2)=2*(b(1)*c(3)-c(1)*b(3));
                 rwx(3)=2*(c(1)*d(4)-d(1)*c(4));        rwy(3)=2*(d(1)*b(4)-b(1)*d(4));          rwz(3)=2*(b(1)*c(4)-c(1)*b(4));
                 rwx(4)=2*(c(2)*d(3)-d(2)*c(3));        rwy(4)=2*(d(2)*b(3)-b(2)*d(3));          rwz(4)=2*(b(2)*c(3)-c(2)*b(3));
                 rwx(5)=2*(c(2)*d(4)-d(2)*c(4));        rwy(5)=2*(d(2)*b(4)-b(2)*d(4));          rwz(5)=2*(b(2)*c(4)-c(2)*b(4));
                 rwx(6)=2*(c(3)*d(4)-d(3)*c(4));        rwy(6)=2*(d(3)*b(4)-b(3)*d(4));          rwz(6)=2*(b(3)*c(4)-c(3)*b(4));
                 %---------------------- Facet Basis Functions ------------
                 wfx(1)=2*zeta(3)*(c(1)*d(2)-c(2)*d(1))+2*zeta(1)*(c(2)*d(3)-c(3)*d(2))+2*zeta(2)*(c(3)*d(1)-c(1)*d(3));    wfy(1)=2*zeta(3)*(d(1)*b(2)-d(2)*b(1))+2*zeta(1)*(d(2)*b(3)-d(3)*b(2))+2*zeta(2)*(d(3)*b(1)-d(1)*b(3));     wfz(1)=2*zeta(3)*(b(1)*c(2)-b(2)*c(1))+2*zeta(1)*(b(2)*c(3)-b(3)*c(2))+2*zeta(2)*(b(3)*c(1)-b(1)*c(3));
                 wfx(2)=2*zeta(3)*(c(2)*d(4)-c(4)*d(2))+2*zeta(2)*(c(4)*d(3)-c(3)*d(4))+2*zeta(4)*(c(3)*d(2)-c(2)*d(3));    wfy(2)=2*zeta(3)*(d(2)*b(4)-d(4)*b(2))+2*zeta(2)*(d(4)*b(3)-d(3)*b(4))+2*zeta(4)*(d(3)*b(2)-d(2)*b(3));     wfz(2)=2*zeta(3)*(b(2)*c(4)-b(4)*c(2))+2*zeta(2)*(b(4)*c(3)-b(3)*c(4))+2*zeta(4)*(b(3)*c(2)-b(2)*c(3));
                 wfx(3)=2*zeta(1)*(c(3)*d(4)-c(4)*d(3))+2*zeta(3)*(c(4)*d(1)-c(1)*d(4))+2*zeta(4)*(c(1)*d(3)-c(3)*d(1));    wfy(3)=2*zeta(1)*(d(3)*b(4)-d(4)*b(3))+2*zeta(3)*(d(4)*b(1)-d(1)*b(4))+2*zeta(4)*(d(1)*b(3)-d(3)*b(1));     wfz(3)=2*zeta(1)*(b(3)*c(4)-b(4)*c(3))+2*zeta(3)*(b(4)*c(1)-b(1)*c(4))+2*zeta(4)*(b(1)*c(3)-b(3)*c(1));
                 wfx(4)=2*zeta(1)*(c(4)*d(2)-c(2)*d(4))+2*zeta(4)*(c(2)*d(1)-c(1)*d(2))+2*zeta(2)*(c(1)*d(4)-c(4)*d(1));    wfy(4)=2*zeta(1)*(d(4)*b(2)-d(2)*b(4))+2*zeta(4)*(d(2)*b(1)-d(1)*b(2))+2*zeta(2)*(d(1)*b(4)-d(4)*b(1));     wfz(4)=2*zeta(1)*(b(4)*c(2)-b(2)*c(4))+2*zeta(4)*(b(2)*c(1)-b(1)*c(2))+2*zeta(2)*(b(1)*c(4)-b(4)*c(1));
                %===================== Basis Function Scaling ============
                 switch TModel.Assembled.E_Scaling
                     case 1,wx=wx.*edgeLengths;     rwx=rwx.*edgeLengths;
                            wy=wy.*edgeLengths;     rwy=rwy.*edgeLengths;
                            wz=wz.*edgeLengths;     rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch TModel.Assembled.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*TModel.Assembled.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*TModel.Assembled.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*TModel.TModel.Assembled.ScalingVec_B(facets(ii).Index);end
                 end
                 %------------------------ (mr^-1)*wf----------------------
                 mwfx(1)=imu(1,1)*wfx(1)+imu(1,2)*wfy(1)+imu(1,3)*wfz(1);   mwfy(1)=imu(2,1)*wfx(1)+imu(2,2)*wfy(1)+imu(2,3)*wfz(1);    mwfz(1)=imu(3,1)*wfx(1)+imu(3,2)*wfy(1)+imu(3,3)*wfz(1);
                 mwfx(2)=imu(1,1)*wfx(2)+imu(1,2)*wfy(2)+imu(1,3)*wfz(2);   mwfy(2)=imu(2,1)*wfx(2)+imu(2,2)*wfy(2)+imu(2,3)*wfz(2);    mwfz(2)=imu(3,1)*wfx(2)+imu(3,2)*wfy(2)+imu(3,3)*wfz(2);
                 mwfx(3)=imu(1,1)*wfx(3)+imu(1,2)*wfy(3)+imu(1,3)*wfz(3);   mwfy(3)=imu(2,1)*wfx(3)+imu(2,2)*wfy(3)+imu(2,3)*wfz(3);    mwfz(3)=imu(3,1)*wfx(3)+imu(3,2)*wfy(3)+imu(3,3)*wfz(3);
                 mwfx(4)=imu(1,1)*wfx(4)+imu(1,2)*wfy(4)+imu(1,3)*wfz(4);   mwfy(4)=imu(2,1)*wfx(4)+imu(2,2)*wfy(4)+imu(2,3)*wfz(4);    mwfz(4)=imu(3,1)*wfx(4)+imu(3,2)*wfy(4)+imu(3,3)*wfz(4);
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
        %======================== Surface Integral Terms ==================
        for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                switch boundary.Type
                    case "ABC"
                        Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(medium.IsDispersive),ZW=medium.WaveImpedance{freqIndex};else,ZW=medium.WaveImpedance;end,Ts=(ZW^-1)*Ts;Ts=Ts*Z0;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "ABB"
                        Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(boundary.Dispersive),beta=boundary.Param(freqIndex);else,beta=boundary.Param;end,Ts=Ts*beta/k0;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "ABZ"
                       if(boundary.Tensor)
                            if(boundary.Dispersive),Z=boundary.Param{freqIndex};
                            else,Z=boundary.Param;
                            end
                            Ts=CalculateTensorIntegral(TModel,Z^-1,element,facets(kk),b,c,d,kk);Ts=Ts*Z0;
                        else
                            if(boundary.Dispersive),Z=boundary.Param(freqIndex);
                            else,Z=boundary.Param;
                            end
                            Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=Ts*Z^-1;Ts=Ts*Z0;
                        end
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "GRA"
                        Tg = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(boundary.Dispersive),cond=boundary.Param(freqIndex);else,cond=boundary.Param;end,Tg=Tg.*cond;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "IBC"
                        if(boundary.Dispersive),cond=boundary.Param(freqIndex);else,cond=boundary.Param;end,Tbc=Tbc.*cond;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "POR"
                        if(boundary.Dispersive)
                            if(boundary.PortParamType==0)
                                Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Ts=Ts*boundary.Param(freqIndex)/k0;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==1),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Zp=boundary.Param(freqIndex);
                                Ts=(Z0/Zp)*Ts;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==2)
                                Zp=boundary.Param{freqIndex};
                                Ts=CalculateTensorIntegral(TModel,Zp^-1,element,facets(kk),b,c,d,kk);
                                Ts=Z0*Ts;Tp=2*1i*Ts;
                            end
                        else
                            if(boundary.PortParamType==0)
                                Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Ts=boundary.Param*Ts/k0;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==1),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Zp=boundary.Param; 
                                Ts=(Z0/Zp)*Ts;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==2),Zp=boundary.Param;Ts=CalculateTensorIntegral(TModel,Zp^-1,element,facets(kk),b,c,d,kk);
                                Ts=Z0*Ts;Tp=2*1i*Ts;
                            end
                        end
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6
                                        for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);
                                                   Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);
                                        end
                                    end
                        end
                end
            end
        end
        %==================== Global Matrix Assembly ======================
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TP(counterA)=dsi*dsj*si*sj*Tp(ii,jj);
                end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;TBB(counterB)=si*sj*Te(ii,jj);end
            end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;AB(counterB)=si*sj*Am(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;FB(counterB)=si*sj*Fm(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj); end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;TBB(counterB)=si*sj*dsi*dsj*Tb(ii,jj);end
            end
        end
    end
end
%--------------------------------------------------------------------------
function [IA,JA,TE,TB,FF,AA,TS,TP,TG,TBC,P,TA,TC,IB,JB,TEB,TBB,FB,AB,PB,TAB,TCB,counterA,counterB] = DirichletExcitation_Bia(varargin),ElectromagneticConstants;GaussianQuadratture3D;
    if(nargin==26)
        IA=varargin{1};         TS=varargin{7};         P=varargin{11};        IB=varargin{14};          PB=varargin{20};           counterA=varargin{23};
        JA=varargin{2};         TP=varargin{8};         TA=varargin{12};       JB=varargin{15};          TAB=varargin{21};          counterB=varargin{24};
        TE=varargin{3};         TG=varargin{9};         TC=varargin{13};       TEB=varargin{16};         TCB=varargin{22};          TModel=varargin{25};
        TB=varargin{4};         TBC=varargin{10};                              TBB=varargin{17};                                    domain=varargin{26};
        FF=varargin{5};                                                        FB=varargin{18};
        AA=varargin{6};                                                        AB=varargin{19};
        medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;ksi=medium.Ksi;zita=medium.Zita;imu=mu^-1;                    
    elseif(nargin==27)
        IA=varargin{1};         TS=varargin{7};         P=varargin{11};        IB=varargin{14};          PB=varargin{20};           counterA=varargin{23};
        JA=varargin{2};         TP=varargin{8};         TA=varargin{12};       JB=varargin{15};          TAB=varargin{21};          counterB=varargin{24};
        TE=varargin{3};         TG=varargin{9};         TC=varargin{13};       TEB=varargin{16};         TCB=varargin{22};          TModel=varargin{25};
        TB=varargin{4};         TBC=varargin{10};                              TBB=varargin{17};                                    domain=varargin{26};
        FF=varargin{5};                                                        FB=varargin{18};                                     freqIndex=varargin{27};
        AA=varargin{6};                                                        AB=varargin{19};
        medium=domain.Medium;
        if(medium.IsDispersive),epsilon=medium.Epsilon{freqIndex};mu=medium.Mu{freqIndex};imu=mu^-1;ksi=medium.Ksi{freqIndex};zita=medium.Zita{freqIndex};
        else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;ksi=medium.Ksi;zita=medium.Zita;
        end
    end,imz=imu*zita;kim=ksi*imu;
    for ie = 1 : numel(domain.Elements),element=TModel.Elements(domain.Elements(ie));
        b=element.Bs;c=element.Cs;d=element.Ds;Ve=element.Volume;edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];
        edgeSigns=[element.EdgeSigns];facetSigns=[element.FacetSigns];edgeLengths=[edges.Length];facetSurfaces=[facets.Surface];
        Tbc=zeros(6,6);Ts=zeros(6,6);Tg=zeros(6,6);Te=zeros(6,6);Tb=zeros(4,4);Fm=zeros(4,6);Am=zeros(6,4);Pm=zeros(6,6);Tam=zeros(6,6);Tcm=zeros(6,4);
        for kt=1:numel(Weights),zeta(1)=Points(1,kt);zeta(2)=Points(2,kt);zeta(3)=Points(3,kt);zeta(4)=Points(4,kt);
                 %---------------------- Edge Basis Functions -------------
                 wx(1)=zeta(1)*b(2)-zeta(2)*b(1);       wy(1)=zeta(1)*c(2)-zeta(2)*c(1);         wz(1)=zeta(1)*d(2)-zeta(2)*d(1);  
                 wx(2)=zeta(1)*b(3)-zeta(3)*b(1);       wy(2)=zeta(1)*c(3)-zeta(3)*c(1);         wz(2)=zeta(1)*d(3)-zeta(3)*d(1);
                 wx(3)=zeta(1)*b(4)-zeta(4)*b(1);       wy(3)=zeta(1)*c(4)-zeta(4)*c(1);         wz(3)=zeta(1)*d(4)-zeta(4)*d(1);
                 wx(4)=zeta(2)*b(3)-zeta(3)*b(2);       wy(4)=zeta(2)*c(3)-zeta(3)*c(2);         wz(4)=zeta(2)*d(3)-zeta(3)*d(2);
                 wx(5)=zeta(2)*b(4)-zeta(4)*b(2);       wy(5)=zeta(2)*c(4)-zeta(4)*c(2);         wz(5)=zeta(2)*d(4)-zeta(4)*d(2);
                 wx(6)=zeta(3)*b(4)-zeta(4)*b(3);       wy(6)=zeta(3)*c(4)-zeta(4)*c(3);         wz(6)=zeta(3)*d(4)-zeta(4)*d(3);
                 %------------------ Edge Basis Functions Rotations -------
                 rwx(1)=2*(c(1)*d(2)-d(1)*c(2));        rwy(1)=2*(d(1)*b(2)-b(1)*d(2));          rwz(1)=2*(b(1)*c(2)-c(1)*b(2));
                 rwx(2)=2*(c(1)*d(3)-d(1)*c(3));        rwy(2)=2*(d(1)*b(3)-b(1)*d(3));          rwz(2)=2*(b(1)*c(3)-c(1)*b(3));
                 rwx(3)=2*(c(1)*d(4)-d(1)*c(4));        rwy(3)=2*(d(1)*b(4)-b(1)*d(4));          rwz(3)=2*(b(1)*c(4)-c(1)*b(4));
                 rwx(4)=2*(c(2)*d(3)-d(2)*c(3));        rwy(4)=2*(d(2)*b(3)-b(2)*d(3));          rwz(4)=2*(b(2)*c(3)-c(2)*b(3));
                 rwx(5)=2*(c(2)*d(4)-d(2)*c(4));        rwy(5)=2*(d(2)*b(4)-b(2)*d(4));          rwz(5)=2*(b(2)*c(4)-c(2)*b(4));
                 rwx(6)=2*(c(3)*d(4)-d(3)*c(4));        rwy(6)=2*(d(3)*b(4)-b(3)*d(4));          rwz(6)=2*(b(3)*c(4)-c(3)*b(4));
                 %---------------------- Facet Basis Functions ------------
                 wfx(1)=2*zeta(3)*(c(1)*d(2)-c(2)*d(1))+2*zeta(1)*(c(2)*d(3)-c(3)*d(2))+2*zeta(2)*(c(3)*d(1)-c(1)*d(3));    wfy(1)=2*zeta(3)*(d(1)*b(2)-d(2)*b(1))+2*zeta(1)*(d(2)*b(3)-d(3)*b(2))+2*zeta(2)*(d(3)*b(1)-d(1)*b(3));     wfz(1)=2*zeta(3)*(b(1)*c(2)-b(2)*c(1))+2*zeta(1)*(b(2)*c(3)-b(3)*c(2))+2*zeta(2)*(b(3)*c(1)-b(1)*c(3));
                 wfx(2)=2*zeta(3)*(c(2)*d(4)-c(4)*d(2))+2*zeta(2)*(c(4)*d(3)-c(3)*d(4))+2*zeta(4)*(c(3)*d(2)-c(2)*d(3));    wfy(2)=2*zeta(3)*(d(2)*b(4)-d(4)*b(2))+2*zeta(2)*(d(4)*b(3)-d(3)*b(4))+2*zeta(4)*(d(3)*b(2)-d(2)*b(3));     wfz(2)=2*zeta(3)*(b(2)*c(4)-b(4)*c(2))+2*zeta(2)*(b(4)*c(3)-b(3)*c(4))+2*zeta(4)*(b(3)*c(2)-b(2)*c(3));
                 wfx(3)=2*zeta(1)*(c(3)*d(4)-c(4)*d(3))+2*zeta(3)*(c(4)*d(1)-c(1)*d(4))+2*zeta(4)*(c(1)*d(3)-c(3)*d(1));    wfy(3)=2*zeta(1)*(d(3)*b(4)-d(4)*b(3))+2*zeta(3)*(d(4)*b(1)-d(1)*b(4))+2*zeta(4)*(d(1)*b(3)-d(3)*b(1));     wfz(3)=2*zeta(1)*(b(3)*c(4)-b(4)*c(3))+2*zeta(3)*(b(4)*c(1)-b(1)*c(4))+2*zeta(4)*(b(1)*c(3)-b(3)*c(1));
                 wfx(4)=2*zeta(1)*(c(4)*d(2)-c(2)*d(4))+2*zeta(4)*(c(2)*d(1)-c(1)*d(2))+2*zeta(2)*(c(1)*d(4)-c(4)*d(1));    wfy(4)=2*zeta(1)*(d(4)*b(2)-d(2)*b(4))+2*zeta(4)*(d(2)*b(1)-d(1)*b(2))+2*zeta(2)*(d(1)*b(4)-d(4)*b(1));     wfz(4)=2*zeta(1)*(b(4)*c(2)-b(2)*c(4))+2*zeta(4)*(b(2)*c(1)-b(1)*c(2))+2*zeta(2)*(b(1)*c(4)-b(4)*c(1));
                 %===================== Basis Function Scaling ============
                 switch TModel.Assembled.E_Scaling
                     case 1,wx=wx.*edgeLengths;     rwx=rwx.*edgeLengths;
                            wy=wy.*edgeLengths;     rwy=rwy.*edgeLengths;
                            wz=wz.*edgeLengths;     rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch TModel.Assembled.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*TModel.Assembled.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*TModel.Assembled.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*TModel.TModel.Assembled.ScalingVec_B(facets(ii).Index);end
                 end
                 %------------------------ (mr^-1)*wf----------------------
                 mwfx(1)=imu(1,1)*wfx(1)+imu(1,2)*wfy(1)+imu(1,3)*wfz(1);   mwfy(1)=imu(2,1)*wfx(1)+imu(2,2)*wfy(1)+imu(2,3)*wfz(1);    mwfz(1)=imu(3,1)*wfx(1)+imu(3,2)*wfy(1)+imu(3,3)*wfz(1);
                 mwfx(2)=imu(1,1)*wfx(2)+imu(1,2)*wfy(2)+imu(1,3)*wfz(2);   mwfy(2)=imu(2,1)*wfx(2)+imu(2,2)*wfy(2)+imu(2,3)*wfz(2);    mwfz(2)=imu(3,1)*wfx(2)+imu(3,2)*wfy(2)+imu(3,3)*wfz(2);
                 mwfx(3)=imu(1,1)*wfx(3)+imu(1,2)*wfy(3)+imu(1,3)*wfz(3);   mwfy(3)=imu(2,1)*wfx(3)+imu(2,2)*wfy(3)+imu(2,3)*wfz(3);    mwfz(3)=imu(3,1)*wfx(3)+imu(3,2)*wfy(3)+imu(3,3)*wfz(3);
                 mwfx(4)=imu(1,1)*wfx(4)+imu(1,2)*wfy(4)+imu(1,3)*wfz(4);   mwfy(4)=imu(2,1)*wfx(4)+imu(2,2)*wfy(4)+imu(2,3)*wfz(4);    mwfz(4)=imu(3,1)*wfx(4)+imu(3,2)*wfy(4)+imu(3,3)*wfz(4);
                %----------- (mr^-1)*zita* Edge Basis Functions -----------
                mzwx(1)=imz(1,1)*wx(1)+imz(1,2)*wy(1)+imz(1,3)*wz(1);       mzwy(1)=imz(2,1)*wx(1)+imz(2,2)*wy(1)+imz(2,3)*wz(1);       mzwz(1)=imz(3,1)*wx(1)+imz(3,2)*wy(1)+imz(3,3)*wz(1);
                mzwx(2)=imz(1,1)*wx(2)+imz(1,2)*wy(2)+imz(1,3)*wz(2);       mzwy(2)=imz(2,1)*wx(2)+imz(2,2)*wy(2)+imz(2,3)*wz(2);       mzwz(2)=imz(3,1)*wx(2)+imz(3,2)*wy(2)+imz(3,3)*wz(2);
                mzwx(3)=imz(1,1)*wx(3)+imz(1,2)*wy(3)+imz(1,3)*wz(3);       mzwy(3)=imz(2,1)*wx(3)+imz(2,2)*wy(3)+imz(2,3)*wz(3);       mzwz(3)=imz(3,1)*wx(3)+imz(3,2)*wy(3)+imz(3,3)*wz(3);
                mzwx(4)=imz(1,1)*wx(4)+imz(1,2)*wy(4)+imz(1,3)*wz(4);       mzwy(4)=imz(2,1)*wx(4)+imz(2,2)*wy(4)+imz(2,3)*wz(4);       mzwz(4)=imz(3,1)*wx(4)+imz(3,2)*wy(4)+imz(3,3)*wz(4);
                mzwx(5)=imz(1,1)*wx(5)+imz(1,2)*wy(5)+imz(1,3)*wz(5);       mzwy(5)=imz(2,1)*wx(5)+imz(2,2)*wy(5)+imz(2,3)*wz(5);       mzwz(5)=imz(3,1)*wx(5)+imz(3,2)*wy(5)+imz(3,3)*wz(5);
                mzwx(6)=imz(1,1)*wx(6)+imz(1,2)*wy(6)+imz(1,3)*wz(6);       mzwy(6)=imz(2,1)*wx(6)+imz(2,2)*wy(6)+imz(2,3)*wz(6);       mzwz(6)=imz(3,1)*wx(6)+imz(3,2)*wy(6)+imz(3,3)*wz(6);
                %------------- ksi *(mr^-1)* FacetBasis Functions ---------
                kmvx(1)=kim(1,1)*wfx(1)+kim(1,2)*wfy(1)+kim(1,3)*wfz(1);    kmvy(1)=kim(2,1)*wfx(1)+kim(2,2)*wfy(1)+kim(2,3)*wfz(1);    kmvz(1)=kim(3,1)*wfx(1)+kim(3,2)*wfy(1)+kim(3,3)*wfz(1);
                kmvx(2)=kim(1,1)*wfx(2)+kim(1,2)*wfy(2)+kim(1,3)*wfz(2);    kmvy(2)=kim(2,1)*wfx(2)+kim(2,2)*wfy(2)+kim(2,3)*wfz(2);    kmvz(2)=kim(3,1)*wfx(2)+kim(3,2)*wfy(2)+kim(3,3)*wfz(2);
                kmvx(3)=kim(1,1)*wfx(3)+kim(1,2)*wfy(3)+kim(1,3)*wfz(3);    kmvy(3)=kim(2,1)*wfx(3)+kim(2,2)*wfy(3)+kim(2,3)*wfz(3);    kmvz(3)=kim(3,1)*wfx(3)+kim(3,2)*wfy(3)+kim(3,3)*wfz(3);
                kmvx(4)=kim(1,1)*wfx(4)+kim(1,2)*wfy(4)+kim(1,3)*wfz(4);    kmvy(4)=kim(2,1)*wfx(4)+kim(2,2)*wfy(4)+kim(2,3)*wfz(4);    kmvz(4)=kim(3,1)*wfx(4)+kim(3,2)*wfy(4)+kim(3,3)*wfz(4);
                %---- ksi*(mr^-1)*zita* Edge Basis Functions --------------
                kmvwx(1)=ksi(1,1)*mzwx(1)+ksi(1,2)*mzwy(1)+ksi(1,3)*mzwz(1);    kmvwy(1)=ksi(2,1)*mzwx(1)+ksi(2,2)*mzwy(1)+ksi(2,3)*mzwz(1);    kmvwz(1)=ksi(3,1)*mzwx(1)+ksi(3,2)*mzwy(1)+ksi(3,3)*mzwz(1);
                kmvwx(2)=ksi(1,1)*mzwx(2)+ksi(1,2)*mzwy(2)+ksi(1,3)*mzwz(2);    kmvwy(2)=ksi(2,1)*mzwx(2)+ksi(2,2)*mzwy(2)+ksi(2,3)*mzwz(2);    kmvwz(2)=ksi(3,1)*mzwx(2)+ksi(3,2)*mzwy(2)+ksi(3,3)*mzwz(2);
                kmvwx(3)=ksi(1,1)*mzwx(3)+ksi(1,2)*mzwy(3)+ksi(1,3)*mzwz(3);    kmvwy(3)=ksi(2,1)*mzwx(3)+ksi(2,2)*mzwy(3)+ksi(2,3)*mzwz(3);    kmvwz(3)=ksi(3,1)*mzwx(3)+ksi(3,2)*mzwy(3)+ksi(3,3)*mzwz(3);
                kmvwx(4)=ksi(1,1)*mzwx(4)+ksi(1,2)*mzwy(4)+ksi(1,3)*mzwz(4);    kmvwy(4)=ksi(2,1)*mzwx(4)+ksi(2,2)*mzwy(4)+ksi(2,3)*mzwz(4);    kmvwz(4)=ksi(3,1)*mzwx(4)+ksi(3,2)*mzwy(4)+ksi(3,3)*mzwz(4);
                kmvwx(5)=ksi(1,1)*mzwx(5)+ksi(1,2)*mzwy(5)+ksi(1,3)*mzwz(5);    kmvwy(5)=ksi(2,1)*mzwx(5)+ksi(2,2)*mzwy(5)+ksi(2,3)*mzwz(5);    kmvwz(5)=ksi(3,1)*mzwx(5)+ksi(3,2)*mzwy(5)+ksi(3,3)*mzwz(5);
                kmvwx(6)=ksi(1,1)*mzwx(6)+ksi(1,2)*mzwy(6)+ksi(1,3)*mzwz(6);    kmvwy(6)=ksi(2,1)*mzwx(6)+ksi(2,2)*mzwy(6)+ksi(2,3)*mzwz(6);    kmvwz(6)=ksi(3,1)*mzwx(6)+ksi(3,2)*mzwy(6)+ksi(3,3)*mzwz(6);
                %=================== Local Matrices =======================
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
         %======================== Surface Integral Terms ==================
        for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                switch boundary.Type
                    case "ABC"
                        Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(medium.IsDispersive),ZW=medium.WaveImpedance{freqIndex};else,ZW=medium.WaveImpedance;end,Ts=(ZW^-1)*Ts;Ts=Ts*Z0;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "ABB"
                        Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(boundary.Dispersive),beta=boundary.Param(freqIndex);else,beta=boundary.Param;end,Ts=Ts*beta/k0;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "ABZ"
                       if(boundary.Tensor)
                            if(boundary.Dispersive),Z=boundary.Param{freqIndex};
                            else,Z=boundary.Param;
                            end
                            Ts=CalculateTensorIntegral(TModel,Z^-1,element,facets(kk),b,c,d,kk);Ts=Ts*Z0;
                        else
                            if(boundary.Dispersive),Z=boundary.Param(freqIndex);
                            else,Z=boundary.Param;
                            end
                            Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=Ts*Z^-1;Ts=Ts*Z0;
                        end
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "GRA"
                        Tg = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(boundary.Dispersive),cond=boundary.Param(freqIndex);else,cond=boundary.Param;end,Tg=Tg.*cond;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "IBC"
                        if(boundary.Dispersive),cond=boundary.Param(freqIndex);else,cond=boundary.Param;end,Tbc=Tbc.*cond;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "POR"
                        if(boundary.Dispersive)
                            if(boundary.PortParamType==0)
                                Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Ts=Ts*boundary.Param(freqIndex)/k0;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==1),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Zp=boundary.Param(freqIndex);
                                Ts=(Z0/Zp)*Ts;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==2)
                                Zp=boundary.Param{freqIndex};
                                Ts=CalculateTensorIntegral(TModel,Zp^-1,element,facets(kk),b,c,d,kk);
                                Ts=Z0*Ts;Tp=2*1i*Ts;
                            end
                        else
                            if(boundary.PortParamType==0)
                                Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Ts=boundary.Param*Ts/k0;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==1),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Zp=boundary.Param; 
                                Ts=(Z0/Zp)*Ts;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==2),Zp=boundary.Param;Ts=CalculateTensorIntegral(TModel,Zp^-1,element,facets(kk),b,c,d,kk);
                                Ts=Z0*Ts;Tp=2*1i*Ts;
                            end
                        end
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6
                                        for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);
                                                   Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);
                                        end
                                    end
                        end
                end
            end
        end
        %==================== Global Matrix Assembly ======================
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TA(counterA)=dsi*dsj*si*sj*Tam(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;P(counterA)=dsi*dsj*si*sj*Pm(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;TEB(counterB)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;TAB(counterB)=dsi*dsj*si*sj*Tam(ii,jj);
                                   counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;PB(counterB)=-dsi*dsj*si*sj*Pm(ii,jj);
                end
            end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TC(counterA)=si*dsi*sj*dsj*Tcm(ii,jj);
                end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterA)=kj;AB(counterB)=si*dsi*sj*dsj*Am(ii,jj);
                                   counterB=counterB+1;IB(counterB)=di;JB(counterA)=kj;TCB(counterB)=si*dsi*sj*dsj*Tcm(ii,jj);

                end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);kj=edges(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IB(counterB)=di;JB(counterB)=kj;FB(counterB)=si*dsi*sj*dsj*Fm(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);kj=facets(jj).KnownIndex;
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj);end
                if(di~=0 && kj~=0),counterB=counterB+1;IA(counterB)=di;JB(counterB)=kj;TBB(counterB)=si*sj*dsi*dsj*Tb(ii,jj);end
            end
        end
     end
end
%========================= Port Excitation Assembly =======================
function [IA,JA,TE,TB,FF,AA,TS,TP,TG,TBC,counterA] =  PortExcitation_Iso(varargin),ElectromagneticConstants;GaussianQuadratture3D;
    if(nargin==13)
        IA=varargin{1};         TS=varargin{7};         counterA=varargin{11};
        JA=varargin{2};         TP=varargin{8};         TModel=varargin{12};
        TE=varargin{3};         TG=varargin{9};         domain=varargin{13};
        TB=varargin{4};         TBC=varargin{10};
        FF=varargin{5};
        AA=varargin{6};
        medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;freq=TModel.Frequency.Frequency;
    elseif(nargin==14)
        IA=varargin{1};         TS=varargin{7};         counterA=varargin{11};
        JA=varargin{2};         TP=varargin{8};         TModel=varargin{12};
        TE=varargin{3};         TG=varargin{9};         domain=varargin{13};
        TB=varargin{4};         TBC=varargin{10};       freqIndex=varargin{14};
        FF=varargin{5};
        AA=varargin{6};
        medium=domain.Medium;freq=TModel.Frequency.Frequency(freqIndex);
        if(medium.IsDispersive),epsilon=medium.Epsilon(freqIndex);mu=medium.Mu(freqIndex);imu=mu^-1;
        else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;
        end
    end,k0=2*pi*freq/c0;
    for ie=1:numel(domain.Elements),element=TModel.Elements(domain.Elements(ie));
        b=element.Bs;c=element.Cs;d=element.Ds;Ve=element.Volume;edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];
        edgeSigns=[element.EdgeSigns];facetSigns=[element.FacetSigns];edgeLengths=[edges.Length];facetSurfaces=[facets.Surface];
        Te=zeros(6,6);Tb=zeros(4,4);Am=zeros(6,4);Fm=zeros(4,6);Tp=zeros(6,6);Ts=zeros(6,6);Tbc=zeros(6,6);Tg=zeros(6,6);
        for kt =1: numel(Weights),zeta(1)=Points(1,kt);zeta(2)=Points(2,kt);zeta(3)=Points(3,kt);zeta(4)=Points(4,kt);
                 %---------------------- Edge Basis Functions -------------
                 wx(1)=zeta(1)*b(2)-zeta(2)*b(1);       wy(1)=zeta(1)*c(2)-zeta(2)*c(1);         wz(1)=zeta(1)*d(2)-zeta(2)*d(1);  
                 wx(2)=zeta(1)*b(3)-zeta(3)*b(1);       wy(2)=zeta(1)*c(3)-zeta(3)*c(1);         wz(2)=zeta(1)*d(3)-zeta(3)*d(1);
                 wx(3)=zeta(1)*b(4)-zeta(4)*b(1);       wy(3)=zeta(1)*c(4)-zeta(4)*c(1);         wz(3)=zeta(1)*d(4)-zeta(4)*d(1);
                 wx(4)=zeta(2)*b(3)-zeta(3)*b(2);       wy(4)=zeta(2)*c(3)-zeta(3)*c(2);         wz(4)=zeta(2)*d(3)-zeta(3)*d(2);
                 wx(5)=zeta(2)*b(4)-zeta(4)*b(2);       wy(5)=zeta(2)*c(4)-zeta(4)*c(2);         wz(5)=zeta(2)*d(4)-zeta(4)*d(2);
                 wx(6)=zeta(3)*b(4)-zeta(4)*b(3);       wy(6)=zeta(3)*c(4)-zeta(4)*c(3);         wz(6)=zeta(3)*d(4)-zeta(4)*d(3);
                 %------------------ Edge Basis Functions Rotations -------
                 rwx(1)=2*(c(1)*d(2)-d(1)*c(2));        rwy(1)=2*(d(1)*b(2)-b(1)*d(2));          rwz(1)=2*(b(1)*c(2)-c(1)*b(2));
                 rwx(2)=2*(c(1)*d(3)-d(1)*c(3));        rwy(2)=2*(d(1)*b(3)-b(1)*d(3));          rwz(2)=2*(b(1)*c(3)-c(1)*b(3));
                 rwx(3)=2*(c(1)*d(4)-d(1)*c(4));        rwy(3)=2*(d(1)*b(4)-b(1)*d(4));          rwz(3)=2*(b(1)*c(4)-c(1)*b(4));
                 rwx(4)=2*(c(2)*d(3)-d(2)*c(3));        rwy(4)=2*(d(2)*b(3)-b(2)*d(3));          rwz(4)=2*(b(2)*c(3)-c(2)*b(3));
                 rwx(5)=2*(c(2)*d(4)-d(2)*c(4));        rwy(5)=2*(d(2)*b(4)-b(2)*d(4));          rwz(5)=2*(b(2)*c(4)-c(2)*b(4));
                 rwx(6)=2*(c(3)*d(4)-d(3)*c(4));        rwy(6)=2*(d(3)*b(4)-b(3)*d(4));          rwz(6)=2*(b(3)*c(4)-c(3)*b(4));
                 %---------------------- Facet Basis Functions ------------
                 wfx(1)=2*zeta(3)*(c(1)*d(2)-c(2)*d(1))+2*zeta(1)*(c(2)*d(3)-c(3)*d(2))+2*zeta(2)*(c(3)*d(1)-c(1)*d(3));    wfy(1)=2*zeta(3)*(d(1)*b(2)-d(2)*b(1))+2*zeta(1)*(d(2)*b(3)-d(3)*b(2))+2*zeta(2)*(d(3)*b(1)-d(1)*b(3));     wfz(1)=2*zeta(3)*(b(1)*c(2)-b(2)*c(1))+2*zeta(1)*(b(2)*c(3)-b(3)*c(2))+2*zeta(2)*(b(3)*c(1)-b(1)*c(3));
                 wfx(2)=2*zeta(3)*(c(2)*d(4)-c(4)*d(2))+2*zeta(2)*(c(4)*d(3)-c(3)*d(4))+2*zeta(4)*(c(3)*d(2)-c(2)*d(3));    wfy(2)=2*zeta(3)*(d(2)*b(4)-d(4)*b(2))+2*zeta(2)*(d(4)*b(3)-d(3)*b(4))+2*zeta(4)*(d(3)*b(2)-d(2)*b(3));     wfz(2)=2*zeta(3)*(b(2)*c(4)-b(4)*c(2))+2*zeta(2)*(b(4)*c(3)-b(3)*c(4))+2*zeta(4)*(b(3)*c(2)-b(2)*c(3));
                 wfx(3)=2*zeta(1)*(c(3)*d(4)-c(4)*d(3))+2*zeta(3)*(c(4)*d(1)-c(1)*d(4))+2*zeta(4)*(c(1)*d(3)-c(3)*d(1));    wfy(3)=2*zeta(1)*(d(3)*b(4)-d(4)*b(3))+2*zeta(3)*(d(4)*b(1)-d(1)*b(4))+2*zeta(4)*(d(1)*b(3)-d(3)*b(1));     wfz(3)=2*zeta(1)*(b(3)*c(4)-b(4)*c(3))+2*zeta(3)*(b(4)*c(1)-b(1)*c(4))+2*zeta(4)*(b(1)*c(3)-b(3)*c(1));
                 wfx(4)=2*zeta(1)*(c(4)*d(2)-c(2)*d(4))+2*zeta(4)*(c(2)*d(1)-c(1)*d(2))+2*zeta(2)*(c(1)*d(4)-c(4)*d(1));    wfy(4)=2*zeta(1)*(d(4)*b(2)-d(2)*b(4))+2*zeta(4)*(d(2)*b(1)-d(1)*b(2))+2*zeta(2)*(d(1)*b(4)-d(4)*b(1));     wfz(4)=2*zeta(1)*(b(4)*c(2)-b(2)*c(4))+2*zeta(4)*(b(2)*c(1)-b(1)*c(2))+2*zeta(2)*(b(1)*c(4)-b(4)*c(1));
                 %===================== Basis Function Scaling ============
                 switch TModel.Assembled.E_Scaling
                     case 1,wx=wx.*edgeLengths;     rwx=rwx.*edgeLengths;
                            wy=wy.*edgeLengths;     rwy=rwy.*edgeLengths;
                            wz=wz.*edgeLengths;     rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch TModel.Assembled.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*TModel.Assembled.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*TModel.Assembled.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*TModel.TModel.Assembled.ScalingVec_B(facets(ii).Index);end
                 end
                 %==================== Local Matrices =====================
                 for ii=1:6,for jj=1:6,Te(ii,jj)=Te(ii,jj)+Weights(kt)*epsilon*(wx(ii)*wx(jj)+wy(ii)*wy(jj)+wz(ii)*wz(jj));end,end
                 for ii=1:4,for jj=1:4,Tb(ii,jj)=Tb(ii,jj)+Weights(kt)*(wfx(ii)*wfx(jj)+wfy(ii)*wfy(jj)+wfz(ii)*wfz(jj));end,end
                 for ii=1:4,for jj=1:6,Fm(ii,jj)=Fm(ii,jj)+Weights(kt)*(wfx(ii)*rwx(jj)+wfy(ii)*rwy(jj)+wfz(ii)*rwz(jj));end,end
                 for ii=1:6,for jj=1:4,Am(ii,jj)=Am(ii,jj)+Weights(kt)*imu*(rwx(ii)*wfx(jj)+rwy(ii)*wfy(jj)+rwz(ii)*wfz(jj));end,end
        end,Te=Te*Ve;Tb=Tb*Ve;Am=Am*Ve;Fm=Fm*Ve;
        %======================== Surface Integral Terms ==================
        for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                switch boundary.Type
                    case "ABC"
                        Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(medium.IsDispersive),ZW=medium.WaveImpedance(freqIndex);else,ZW=medium.WaveImpedance;end,Ts=(ZW^-1)*Ts;Ts=Ts*Z0;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "ABB"
                        Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(boundary.Dispersive),beta=boundary.Param(freqIndex);else,beta=boundary.Param;end,Ts=Ts*beta/k0;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "ABZ"
                       if(boundary.Tensor)
                            if(boundary.Dispersive),Z=boundary.Param{freqIndex};
                            else,Z=boundary.Param;
                            end
                            Ts=CalculateTensorIntegral(TModel,Z^-1,element,facets(kk),b,c,d,kk);Ts=Ts*Z0;
                        else
                            if(boundary.Dispersive),Z=boundary.Param(freqIndex);
                            else,Z=boundary.Param;
                            end
                            Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=Ts*Z^-1;Ts=Ts*Z0;
                        end
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "GRA"
                        Tg = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(boundary.Dispersive),cond=boundary.Param(freqIndex);else,cond=boundary.Param;end,Tg=Tg.*cond;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "IBC"
                        if(boundary.Dispersive),cond=boundary.Param(freqIndex);else,cond=boundary.Param;end,Tbc=Tbc.*cond;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "POR"
                        if(boundary.Dispersive)
                            if(boundary.PortParamType==0)
                                Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Ts=Ts*boundary.Param(freqIndex)/k0;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==1),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Zp=boundary.Param(freqIndex);
                                Ts=(Z0/Zp)*Ts;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==2)
                                Zp=boundary.Param{freqIndex};
                                Ts=CalculateTensorIntegral(TModel,Zp^-1,element,facets(kk),b,c,d,kk);
                                Ts=Z0*Ts;Tp=2*1i*Ts;
                            end
                        else
                            if(boundary.PortParamType==0)
                                Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Ts=boundary.Param*Ts/k0;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==1),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Zp=boundary.Param; 
                                Ts=(Z0/Zp)*Ts;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==2),Zp=boundary.Param;Ts=CalculateTensorIntegral(TModel,Zp^-1,element,facets(kk),b,c,d,kk);
                                Ts=Z0*Ts;Tp=2*1i*Ts;
                            end
                        end
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6
                                        for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);
                                                   Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);
                                        end
                                    end
                        end
                end
            end
        end
        %==================== Global Matrix Assembly ======================
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TP(counterA)=dsi*dsj*si*sj*Tp(ii,jj);
                end
            end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj); end
            end
        end
    end
end
%--------------------------------------------------------------------------
function [IA,JA,TE,TB,FF,AA,TS,TP,TG,TBC,counterA] =  PortExcitation_Ani(varargin),ElectromagneticConstants;GaussianQuadratture3D;
    if(nargin==13)
        IA=varargin{1};         TS=varargin{7};         counterA=varargin{11};
        JA=varargin{2};         TP=varargin{8};         TModel=varargin{12};
        TE=varargin{3};         TG=varargin{9};         domain=varargin{13};
        TB=varargin{4};         TBC=varargin{10};
        FF=varargin{5};
        AA=varargin{6};
        medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;freq=TModel.Frequency.Frequency;
    elseif(nargin==14)
        IA=varargin{1};         TS=varargin{7};         counterA=varargin{11};
        JA=varargin{2};         TP=varargin{8};         TModel=varargin{12};
        TE=varargin{3};         TG=varargin{9};         domain=varargin{13};
        TB=varargin{4};         TBC=varargin{10};       freqIndex=varargin{14};
        FF=varargin{5};
        AA=varargin{6};
        medium=domain.Medium;freq=TModel.Frequency.Frequency(freqIndex);
        if(medium.IsDispersive),epsilon=medium.Epsilon{freqIndex};mu=medium.Mu{freqIndex};imu=mu^-1;
        else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;
        end
    end,k0=2*pi*freq/c0;
    for ie=1:numel(domain.Elements),element=TModel.Elements(domain.Elements(ie));
        b=element.Bs;c=element.Cs;d=element.Ds;Ve=element.Volume;edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];
        edgeSigns=[element.EdgeSigns];facetSigns=[element.FacetSigns];edgeLengths=[edges.Length];facetSurfaces=[facets.Surface];
        Te=zeros(6,6);Tb=zeros(4,4);Am=zeros(6,4);Fm=zeros(4,6);Tp=zeros(6,6);Ts=zeros(6,6);Tbc=zeros(6,6);Tg=zeros(6,6);
        for kt =1: numel(Weights),zeta(1)=Points(1,kt);zeta(2)=Points(2,kt);zeta(3)=Points(3,kt);zeta(4)=Points(4,kt);
                 %---------------------- Edge Basis Functions -------------
                 wx(1)=zeta(1)*b(2)-zeta(2)*b(1);       wy(1)=zeta(1)*c(2)-zeta(2)*c(1);         wz(1)=zeta(1)*d(2)-zeta(2)*d(1);  
                 wx(2)=zeta(1)*b(3)-zeta(3)*b(1);       wy(2)=zeta(1)*c(3)-zeta(3)*c(1);         wz(2)=zeta(1)*d(3)-zeta(3)*d(1);
                 wx(3)=zeta(1)*b(4)-zeta(4)*b(1);       wy(3)=zeta(1)*c(4)-zeta(4)*c(1);         wz(3)=zeta(1)*d(4)-zeta(4)*d(1);
                 wx(4)=zeta(2)*b(3)-zeta(3)*b(2);       wy(4)=zeta(2)*c(3)-zeta(3)*c(2);         wz(4)=zeta(2)*d(3)-zeta(3)*d(2);
                 wx(5)=zeta(2)*b(4)-zeta(4)*b(2);       wy(5)=zeta(2)*c(4)-zeta(4)*c(2);         wz(5)=zeta(2)*d(4)-zeta(4)*d(2);
                 wx(6)=zeta(3)*b(4)-zeta(4)*b(3);       wy(6)=zeta(3)*c(4)-zeta(4)*c(3);         wz(6)=zeta(3)*d(4)-zeta(4)*d(3);
                 %------------------ Edge Basis Functions Rotations -------
                 rwx(1)=2*(c(1)*d(2)-d(1)*c(2));        rwy(1)=2*(d(1)*b(2)-b(1)*d(2));          rwz(1)=2*(b(1)*c(2)-c(1)*b(2));
                 rwx(2)=2*(c(1)*d(3)-d(1)*c(3));        rwy(2)=2*(d(1)*b(3)-b(1)*d(3));          rwz(2)=2*(b(1)*c(3)-c(1)*b(3));
                 rwx(3)=2*(c(1)*d(4)-d(1)*c(4));        rwy(3)=2*(d(1)*b(4)-b(1)*d(4));          rwz(3)=2*(b(1)*c(4)-c(1)*b(4));
                 rwx(4)=2*(c(2)*d(3)-d(2)*c(3));        rwy(4)=2*(d(2)*b(3)-b(2)*d(3));          rwz(4)=2*(b(2)*c(3)-c(2)*b(3));
                 rwx(5)=2*(c(2)*d(4)-d(2)*c(4));        rwy(5)=2*(d(2)*b(4)-b(2)*d(4));          rwz(5)=2*(b(2)*c(4)-c(2)*b(4));
                 rwx(6)=2*(c(3)*d(4)-d(3)*c(4));        rwy(6)=2*(d(3)*b(4)-b(3)*d(4));          rwz(6)=2*(b(3)*c(4)-c(3)*b(4));
                 %---------------------- Facet Basis Functions ------------
                 wfx(1)=2*zeta(3)*(c(1)*d(2)-c(2)*d(1))+2*zeta(1)*(c(2)*d(3)-c(3)*d(2))+2*zeta(2)*(c(3)*d(1)-c(1)*d(3));    wfy(1)=2*zeta(3)*(d(1)*b(2)-d(2)*b(1))+2*zeta(1)*(d(2)*b(3)-d(3)*b(2))+2*zeta(2)*(d(3)*b(1)-d(1)*b(3));     wfz(1)=2*zeta(3)*(b(1)*c(2)-b(2)*c(1))+2*zeta(1)*(b(2)*c(3)-b(3)*c(2))+2*zeta(2)*(b(3)*c(1)-b(1)*c(3));
                 wfx(2)=2*zeta(3)*(c(2)*d(4)-c(4)*d(2))+2*zeta(2)*(c(4)*d(3)-c(3)*d(4))+2*zeta(4)*(c(3)*d(2)-c(2)*d(3));    wfy(2)=2*zeta(3)*(d(2)*b(4)-d(4)*b(2))+2*zeta(2)*(d(4)*b(3)-d(3)*b(4))+2*zeta(4)*(d(3)*b(2)-d(2)*b(3));     wfz(2)=2*zeta(3)*(b(2)*c(4)-b(4)*c(2))+2*zeta(2)*(b(4)*c(3)-b(3)*c(4))+2*zeta(4)*(b(3)*c(2)-b(2)*c(3));
                 wfx(3)=2*zeta(1)*(c(3)*d(4)-c(4)*d(3))+2*zeta(3)*(c(4)*d(1)-c(1)*d(4))+2*zeta(4)*(c(1)*d(3)-c(3)*d(1));    wfy(3)=2*zeta(1)*(d(3)*b(4)-d(4)*b(3))+2*zeta(3)*(d(4)*b(1)-d(1)*b(4))+2*zeta(4)*(d(1)*b(3)-d(3)*b(1));     wfz(3)=2*zeta(1)*(b(3)*c(4)-b(4)*c(3))+2*zeta(3)*(b(4)*c(1)-b(1)*c(4))+2*zeta(4)*(b(1)*c(3)-b(3)*c(1));
                 wfx(4)=2*zeta(1)*(c(4)*d(2)-c(2)*d(4))+2*zeta(4)*(c(2)*d(1)-c(1)*d(2))+2*zeta(2)*(c(1)*d(4)-c(4)*d(1));    wfy(4)=2*zeta(1)*(d(4)*b(2)-d(2)*b(4))+2*zeta(4)*(d(2)*b(1)-d(1)*b(2))+2*zeta(2)*(d(1)*b(4)-d(4)*b(1));     wfz(4)=2*zeta(1)*(b(4)*c(2)-b(2)*c(4))+2*zeta(4)*(b(2)*c(1)-b(1)*c(2))+2*zeta(2)*(b(1)*c(4)-b(4)*c(1));
                 %===================== Basis Function Scaling ============
                 switch TModel.Assembled.E_Scaling
                     case 1,wx=wx.*edgeLengths;     rwx=rwx.*edgeLengths;
                            wy=wy.*edgeLengths;     rwy=rwy.*edgeLengths;
                            wz=wz.*edgeLengths;     rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch TModel.Assembled.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*TModel.Assembled.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*TModel.Assembled.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*TModel.TModel.Assembled.ScalingVec_B(facets(ii).Index);end
                 end
                 %------------------------ (mr^-1)*wf----------------------
                 mwfx(1)=imu(1,1)*wfx(1)+imu(1,2)*wfy(1)+imu(1,3)*wfz(1);   mwfy(1)=imu(2,1)*wfx(1)+imu(2,2)*wfy(1)+imu(2,3)*wfz(1);    mwfz(1)=imu(3,1)*wfx(1)+imu(3,2)*wfy(1)+imu(3,3)*wfz(1);
                 mwfx(2)=imu(1,1)*wfx(2)+imu(1,2)*wfy(2)+imu(1,3)*wfz(2);   mwfy(2)=imu(2,1)*wfx(2)+imu(2,2)*wfy(2)+imu(2,3)*wfz(2);    mwfz(2)=imu(3,1)*wfx(2)+imu(3,2)*wfy(2)+imu(3,3)*wfz(2);
                 mwfx(3)=imu(1,1)*wfx(3)+imu(1,2)*wfy(3)+imu(1,3)*wfz(3);   mwfy(3)=imu(2,1)*wfx(3)+imu(2,2)*wfy(3)+imu(2,3)*wfz(3);    mwfz(3)=imu(3,1)*wfx(3)+imu(3,2)*wfy(3)+imu(3,3)*wfz(3);
                 mwfx(4)=imu(1,1)*wfx(4)+imu(1,2)*wfy(4)+imu(1,3)*wfz(4);   mwfy(4)=imu(2,1)*wfx(4)+imu(2,2)*wfy(4)+imu(2,3)*wfz(4);    mwfz(4)=imu(3,1)*wfx(4)+imu(3,2)*wfy(4)+imu(3,3)*wfz(4);
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
        %======================== Surface Integral Terms ==================
        for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                switch boundary.Type
                    case "ABC"
                        Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(medium.IsDispersive),ZW=medium.WaveImpedance{freqIndex};else,ZW=medium.WaveImpedance;end,Ts=(ZW^-1)*Ts;Ts=Ts*Z0;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "ABB"
                        Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(boundary.Dispersive),beta=boundary.Param(freqIndex);else,beta=boundary.Param;end,Ts=Ts*beta/k0;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "ABZ"
                       if(boundary.Tensor)
                            if(boundary.Dispersive),Z=boundary.Param{freqIndex};
                            else,Z=boundary.Param;
                            end
                            Ts=CalculateTensorIntegral(TModel,Z^-1,element,facets(kk),b,c,d,kk);Ts=Ts*Z0;
                        else
                            if(boundary.Dispersive),Z=boundary.Param(freqIndex);
                            else,Z=boundary.Param;
                            end
                            Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=Ts*Z^-1;Ts=Ts*Z0;
                        end
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "GRA"
                        Tg = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(boundary.Dispersive),cond=boundary.Param(freqIndex);else,cond=boundary.Param;end,Tg=Tg.*cond;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "IBC"
                        if(boundary.Dispersive),cond=boundary.Param(freqIndex);else,cond=boundary.Param;end,Tbc=Tbc.*cond;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "POR"
                        if(boundary.Dispersive)
                            if(boundary.PortParamType==0)
                                Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Ts=Ts*boundary.Param(freqIndex)/k0;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==1),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Zp=boundary.Param(freqIndex);
                                Ts=(Z0/Zp)*Ts;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==2)
                                Zp=boundary.Param{freqIndex};
                                Ts=CalculateTensorIntegral(TModel,Zp^-1,element,facets(kk),b,c,d,kk);
                                Ts=Z0*Ts;Tp=2*1i*Ts;
                            end
                        else
                            if(boundary.PortParamType==0)
                                Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Ts=boundary.Param*Ts/k0;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==1),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Zp=boundary.Param; 
                                Ts=(Z0/Zp)*Ts;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==2),Zp=boundary.Param;Ts=CalculateTensorIntegral(TModel,Zp^-1,element,facets(kk),b,c,d,kk);
                                Ts=Z0*Ts;Tp=2*1i*Ts;
                            end
                        end
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6
                                        for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);
                                                   Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);
                                        end
                                    end
                        end
                end
            end
        end
        %==================== Global Matrix Assembly ======================
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TE(counterA)=dsi*dsj*si*sj*Te(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TS(counterA)=dsi*dsj*si*sj*Ts(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TG(counterA)=dsi*dsj*si*sj*Tg(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TBC(counterA)=dsi*dsj*si*sj*Tbc(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TP(counterA)=dsi*dsj*si*sj*Tp(ii,jj);
                end
            end
        end
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);end
             end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj); end
            end
        end
    end
end
%--------------------------------------------------------------------------
function [IA,JA,TE,TB,FF,AA,TS,TP,TG,TBC,P,TA,TC,counterA] = PortExcitation_Bia(varargin),ElectromagneticConstants;GaussianQuadratture3D;
    if(nargin==16)
        IA=varargin{1};         TS=varargin{7};         P=varargin{11};          counterA=varargin{14};
        JA=varargin{2};         TP=varargin{8};         TA=varargin{12};         TModel=varargin{15};     
        TE=varargin{3};         TG=varargin{9};         TC=varargin{13};         domain=varargin{16};        
        TB=varargin{4};         TBC=varargin{10};                                                                 
        FF=varargin{5};                                                    
        AA=varargin{6};                                                        
        medium=domain.Medium;epsilon=medium.Epsilon;mu=medium.Mu;ksi=medium.Ksi;zita=medium.Zita;imu=mu^-1;                    
    elseif(nargin==17)
        IA=varargin{1};         TS=varargin{7};         P=varargin{11};          counterA=varargin{14};
        JA=varargin{2};         TP=varargin{8};         TA=varargin{12};         TModel=varargin{15};     
        TE=varargin{3};         TG=varargin{9};         TC=varargin{13};         domain=varargin{16};        
        TB=varargin{4};         TBC=varargin{10};                                freqIndex=varargin{17};                                                         
        FF=varargin{5};                                                    
        AA=varargin{6};        
        medium=domain.Medium;
        if(medium.IsDispersive),epsilon=medium.Epsilon{freqIndex};mu=medium.Mu{freqIndex};imu=mu^-1;ksi=medium.Ksi{freqIndex};zita=medium.Zita{freqIndex};
        else,epsilon=medium.Epsilon;mu=medium.Mu;imu=mu^-1;ksi=medium.Ksi;zita=medium.Zita;
        end
    end,imz=imu*zita;kim=ksi*imu;
    for ie = 1 : numel(domain.Elements),element=TModel.Elements(domain.Elements(ie));
        b=element.Bs;c=element.Cs;d=element.Ds;Ve=element.Volume;edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];
        edgeSigns=[element.EdgeSigns];facetSigns=[element.FacetSigns];edgeLengths=[edges.Length];facetSurfaces=[facets.Surface];
        Tbc=zeros(6,6);Ts=zeros(6,6);Tg=zeros(6,6);Te=zeros(6,6);Tb=zeros(4,4);Fm=zeros(4,6);Am=zeros(6,4);Pm=zeros(6,6);Tam=zeros(6,6);Tcm=zeros(6,4);
        for kt=1:numel(Weights),zeta(1)=Points(1,kt);zeta(2)=Points(2,kt);zeta(3)=Points(3,kt);zeta(4)=Points(4,kt);
                 %---------------------- Edge Basis Functions -------------
                 wx(1)=zeta(1)*b(2)-zeta(2)*b(1);       wy(1)=zeta(1)*c(2)-zeta(2)*c(1);         wz(1)=zeta(1)*d(2)-zeta(2)*d(1);  
                 wx(2)=zeta(1)*b(3)-zeta(3)*b(1);       wy(2)=zeta(1)*c(3)-zeta(3)*c(1);         wz(2)=zeta(1)*d(3)-zeta(3)*d(1);
                 wx(3)=zeta(1)*b(4)-zeta(4)*b(1);       wy(3)=zeta(1)*c(4)-zeta(4)*c(1);         wz(3)=zeta(1)*d(4)-zeta(4)*d(1);
                 wx(4)=zeta(2)*b(3)-zeta(3)*b(2);       wy(4)=zeta(2)*c(3)-zeta(3)*c(2);         wz(4)=zeta(2)*d(3)-zeta(3)*d(2);
                 wx(5)=zeta(2)*b(4)-zeta(4)*b(2);       wy(5)=zeta(2)*c(4)-zeta(4)*c(2);         wz(5)=zeta(2)*d(4)-zeta(4)*d(2);
                 wx(6)=zeta(3)*b(4)-zeta(4)*b(3);       wy(6)=zeta(3)*c(4)-zeta(4)*c(3);         wz(6)=zeta(3)*d(4)-zeta(4)*d(3);
                 %------------------ Edge Basis Functions Rotations -------
                 rwx(1)=2*(c(1)*d(2)-d(1)*c(2));        rwy(1)=2*(d(1)*b(2)-b(1)*d(2));          rwz(1)=2*(b(1)*c(2)-c(1)*b(2));
                 rwx(2)=2*(c(1)*d(3)-d(1)*c(3));        rwy(2)=2*(d(1)*b(3)-b(1)*d(3));          rwz(2)=2*(b(1)*c(3)-c(1)*b(3));
                 rwx(3)=2*(c(1)*d(4)-d(1)*c(4));        rwy(3)=2*(d(1)*b(4)-b(1)*d(4));          rwz(3)=2*(b(1)*c(4)-c(1)*b(4));
                 rwx(4)=2*(c(2)*d(3)-d(2)*c(3));        rwy(4)=2*(d(2)*b(3)-b(2)*d(3));          rwz(4)=2*(b(2)*c(3)-c(2)*b(3));
                 rwx(5)=2*(c(2)*d(4)-d(2)*c(4));        rwy(5)=2*(d(2)*b(4)-b(2)*d(4));          rwz(5)=2*(b(2)*c(4)-c(2)*b(4));
                 rwx(6)=2*(c(3)*d(4)-d(3)*c(4));        rwy(6)=2*(d(3)*b(4)-b(3)*d(4));          rwz(6)=2*(b(3)*c(4)-c(3)*b(4));
                 %---------------------- Facet Basis Functions ------------
                 wfx(1)=2*zeta(3)*(c(1)*d(2)-c(2)*d(1))+2*zeta(1)*(c(2)*d(3)-c(3)*d(2))+2*zeta(2)*(c(3)*d(1)-c(1)*d(3));    wfy(1)=2*zeta(3)*(d(1)*b(2)-d(2)*b(1))+2*zeta(1)*(d(2)*b(3)-d(3)*b(2))+2*zeta(2)*(d(3)*b(1)-d(1)*b(3));     wfz(1)=2*zeta(3)*(b(1)*c(2)-b(2)*c(1))+2*zeta(1)*(b(2)*c(3)-b(3)*c(2))+2*zeta(2)*(b(3)*c(1)-b(1)*c(3));
                 wfx(2)=2*zeta(3)*(c(2)*d(4)-c(4)*d(2))+2*zeta(2)*(c(4)*d(3)-c(3)*d(4))+2*zeta(4)*(c(3)*d(2)-c(2)*d(3));    wfy(2)=2*zeta(3)*(d(2)*b(4)-d(4)*b(2))+2*zeta(2)*(d(4)*b(3)-d(3)*b(4))+2*zeta(4)*(d(3)*b(2)-d(2)*b(3));     wfz(2)=2*zeta(3)*(b(2)*c(4)-b(4)*c(2))+2*zeta(2)*(b(4)*c(3)-b(3)*c(4))+2*zeta(4)*(b(3)*c(2)-b(2)*c(3));
                 wfx(3)=2*zeta(1)*(c(3)*d(4)-c(4)*d(3))+2*zeta(3)*(c(4)*d(1)-c(1)*d(4))+2*zeta(4)*(c(1)*d(3)-c(3)*d(1));    wfy(3)=2*zeta(1)*(d(3)*b(4)-d(4)*b(3))+2*zeta(3)*(d(4)*b(1)-d(1)*b(4))+2*zeta(4)*(d(1)*b(3)-d(3)*b(1));     wfz(3)=2*zeta(1)*(b(3)*c(4)-b(4)*c(3))+2*zeta(3)*(b(4)*c(1)-b(1)*c(4))+2*zeta(4)*(b(1)*c(3)-b(3)*c(1));
                 wfx(4)=2*zeta(1)*(c(4)*d(2)-c(2)*d(4))+2*zeta(4)*(c(2)*d(1)-c(1)*d(2))+2*zeta(2)*(c(1)*d(4)-c(4)*d(1));    wfy(4)=2*zeta(1)*(d(4)*b(2)-d(2)*b(4))+2*zeta(4)*(d(2)*b(1)-d(1)*b(2))+2*zeta(2)*(d(1)*b(4)-d(4)*b(1));     wfz(4)=2*zeta(1)*(b(4)*c(2)-b(2)*c(4))+2*zeta(4)*(b(2)*c(1)-b(1)*c(2))+2*zeta(2)*(b(1)*c(4)-b(4)*c(1));
                 %===================== Basis Function Scaling ============
                 switch TModel.Assembled.E_Scaling
                     case 1,wx=wx.*edgeLengths;     rwx=rwx.*edgeLengths;
                            wy=wy.*edgeLengths;     rwy=rwy.*edgeLengths;
                            wz=wz.*edgeLengths;     rwz=rwz.*edgeLengths;
                     case 2,for ii=1:6,wx(ii)=wx(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);wy(ii)=wy(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);wz(ii)=wz(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);
                                rwx(ii)=rwx(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);rwy(ii)=rwy(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);rwz(ii)=rwz(ii)*TModel.Assembled.ScalingVec_E(edges(ii).Index);
                            end
                 end
                 switch TModel.Assembled.B_Scaling
                     case 1,wfx=wfx.*facetSurfaces;wfy=wfy.*facetSurfaces;wfz=wfz.*facetSurfaces;
                     case 2,for ii=1:4,wfx(ii)=wfx(ii)*TModel.Assembled.ScalingVec_B(facets(ii).Index);wfy(ii)=wfy(ii)*TModel.Assembled.ScalingVec_B(facets(ii).Index);wfz(ii)=wfz(ii)*TModel.TModel.Assembled.ScalingVec_B(facets(ii).Index);end
                 end
                 %------------------------ (mr^-1)*wf----------------------
                 mwfx(1)=imu(1,1)*wfx(1)+imu(1,2)*wfy(1)+imu(1,3)*wfz(1);   mwfy(1)=imu(2,1)*wfx(1)+imu(2,2)*wfy(1)+imu(2,3)*wfz(1);    mwfz(1)=imu(3,1)*wfx(1)+imu(3,2)*wfy(1)+imu(3,3)*wfz(1);
                 mwfx(2)=imu(1,1)*wfx(2)+imu(1,2)*wfy(2)+imu(1,3)*wfz(2);   mwfy(2)=imu(2,1)*wfx(2)+imu(2,2)*wfy(2)+imu(2,3)*wfz(2);    mwfz(2)=imu(3,1)*wfx(2)+imu(3,2)*wfy(2)+imu(3,3)*wfz(2);
                 mwfx(3)=imu(1,1)*wfx(3)+imu(1,2)*wfy(3)+imu(1,3)*wfz(3);   mwfy(3)=imu(2,1)*wfx(3)+imu(2,2)*wfy(3)+imu(2,3)*wfz(3);    mwfz(3)=imu(3,1)*wfx(3)+imu(3,2)*wfy(3)+imu(3,3)*wfz(3);
                 mwfx(4)=imu(1,1)*wfx(4)+imu(1,2)*wfy(4)+imu(1,3)*wfz(4);   mwfy(4)=imu(2,1)*wfx(4)+imu(2,2)*wfy(4)+imu(2,3)*wfz(4);    mwfz(4)=imu(3,1)*wfx(4)+imu(3,2)*wfy(4)+imu(3,3)*wfz(4);
                %----------- (mr^-1)*zita* Edge Basis Functions -----------
                mzwx(1)=imz(1,1)*wx(1)+imz(1,2)*wy(1)+imz(1,3)*wz(1);       mzwy(1)=imz(2,1)*wx(1)+imz(2,2)*wy(1)+imz(2,3)*wz(1);       mzwz(1)=imz(3,1)*wx(1)+imz(3,2)*wy(1)+imz(3,3)*wz(1);
                mzwx(2)=imz(1,1)*wx(2)+imz(1,2)*wy(2)+imz(1,3)*wz(2);       mzwy(2)=imz(2,1)*wx(2)+imz(2,2)*wy(2)+imz(2,3)*wz(2);       mzwz(2)=imz(3,1)*wx(2)+imz(3,2)*wy(2)+imz(3,3)*wz(2);
                mzwx(3)=imz(1,1)*wx(3)+imz(1,2)*wy(3)+imz(1,3)*wz(3);       mzwy(3)=imz(2,1)*wx(3)+imz(2,2)*wy(3)+imz(2,3)*wz(3);       mzwz(3)=imz(3,1)*wx(3)+imz(3,2)*wy(3)+imz(3,3)*wz(3);
                mzwx(4)=imz(1,1)*wx(4)+imz(1,2)*wy(4)+imz(1,3)*wz(4);       mzwy(4)=imz(2,1)*wx(4)+imz(2,2)*wy(4)+imz(2,3)*wz(4);       mzwz(4)=imz(3,1)*wx(4)+imz(3,2)*wy(4)+imz(3,3)*wz(4);
                mzwx(5)=imz(1,1)*wx(5)+imz(1,2)*wy(5)+imz(1,3)*wz(5);       mzwy(5)=imz(2,1)*wx(5)+imz(2,2)*wy(5)+imz(2,3)*wz(5);       mzwz(5)=imz(3,1)*wx(5)+imz(3,2)*wy(5)+imz(3,3)*wz(5);
                mzwx(6)=imz(1,1)*wx(6)+imz(1,2)*wy(6)+imz(1,3)*wz(6);       mzwy(6)=imz(2,1)*wx(6)+imz(2,2)*wy(6)+imz(2,3)*wz(6);       mzwz(6)=imz(3,1)*wx(6)+imz(3,2)*wy(6)+imz(3,3)*wz(6);
                %------------- ksi *(mr^-1)* FacetBasis Functions ---------
                kmvx(1)=kim(1,1)*wfx(1)+kim(1,2)*wfy(1)+kim(1,3)*wfz(1);    kmvy(1)=kim(2,1)*wfx(1)+kim(2,2)*wfy(1)+kim(2,3)*wfz(1);    kmvz(1)=kim(3,1)*wfx(1)+kim(3,2)*wfy(1)+kim(3,3)*wfz(1);
                kmvx(2)=kim(1,1)*wfx(2)+kim(1,2)*wfy(2)+kim(1,3)*wfz(2);    kmvy(2)=kim(2,1)*wfx(2)+kim(2,2)*wfy(2)+kim(2,3)*wfz(2);    kmvz(2)=kim(3,1)*wfx(2)+kim(3,2)*wfy(2)+kim(3,3)*wfz(2);
                kmvx(3)=kim(1,1)*wfx(3)+kim(1,2)*wfy(3)+kim(1,3)*wfz(3);    kmvy(3)=kim(2,1)*wfx(3)+kim(2,2)*wfy(3)+kim(2,3)*wfz(3);    kmvz(3)=kim(3,1)*wfx(3)+kim(3,2)*wfy(3)+kim(3,3)*wfz(3);
                kmvx(4)=kim(1,1)*wfx(4)+kim(1,2)*wfy(4)+kim(1,3)*wfz(4);    kmvy(4)=kim(2,1)*wfx(4)+kim(2,2)*wfy(4)+kim(2,3)*wfz(4);    kmvz(4)=kim(3,1)*wfx(4)+kim(3,2)*wfy(4)+kim(3,3)*wfz(4);
                %---- ksi*(mr^-1)*zita* Edge Basis Functions --------------
                kmvwx(1)=ksi(1,1)*mzwx(1)+ksi(1,2)*mzwy(1)+ksi(1,3)*mzwz(1);    kmvwy(1)=ksi(2,1)*mzwx(1)+ksi(2,2)*mzwy(1)+ksi(2,3)*mzwz(1);    kmvwz(1)=ksi(3,1)*mzwx(1)+ksi(3,2)*mzwy(1)+ksi(3,3)*mzwz(1);
                kmvwx(2)=ksi(1,1)*mzwx(2)+ksi(1,2)*mzwy(2)+ksi(1,3)*mzwz(2);    kmvwy(2)=ksi(2,1)*mzwx(2)+ksi(2,2)*mzwy(2)+ksi(2,3)*mzwz(2);    kmvwz(2)=ksi(3,1)*mzwx(2)+ksi(3,2)*mzwy(2)+ksi(3,3)*mzwz(2);
                kmvwx(3)=ksi(1,1)*mzwx(3)+ksi(1,2)*mzwy(3)+ksi(1,3)*mzwz(3);    kmvwy(3)=ksi(2,1)*mzwx(3)+ksi(2,2)*mzwy(3)+ksi(2,3)*mzwz(3);    kmvwz(3)=ksi(3,1)*mzwx(3)+ksi(3,2)*mzwy(3)+ksi(3,3)*mzwz(3);
                kmvwx(4)=ksi(1,1)*mzwx(4)+ksi(1,2)*mzwy(4)+ksi(1,3)*mzwz(4);    kmvwy(4)=ksi(2,1)*mzwx(4)+ksi(2,2)*mzwy(4)+ksi(2,3)*mzwz(4);    kmvwz(4)=ksi(3,1)*mzwx(4)+ksi(3,2)*mzwy(4)+ksi(3,3)*mzwz(4);
                kmvwx(5)=ksi(1,1)*mzwx(5)+ksi(1,2)*mzwy(5)+ksi(1,3)*mzwz(5);    kmvwy(5)=ksi(2,1)*mzwx(5)+ksi(2,2)*mzwy(5)+ksi(2,3)*mzwz(5);    kmvwz(5)=ksi(3,1)*mzwx(5)+ksi(3,2)*mzwy(5)+ksi(3,3)*mzwz(5);
                kmvwx(6)=ksi(1,1)*mzwx(6)+ksi(1,2)*mzwy(6)+ksi(1,3)*mzwz(6);    kmvwy(6)=ksi(2,1)*mzwx(6)+ksi(2,2)*mzwy(6)+ksi(2,3)*mzwz(6);    kmvwz(6)=ksi(3,1)*mzwx(6)+ksi(3,2)*mzwy(6)+ksi(3,3)*mzwz(6);
                %=================== Local Matrices =======================
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
         %======================== Surface Integral Terms ==================
        for kk=1:4
            if(~isempty(facets(kk).OnBoundary)),boundary=TModel.Boundaries(facets(kk).OnBoundary);
                switch boundary.Type
                    case "ABC"
                        Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(medium.IsDispersive),ZW=medium.WaveImpedance{freqIndex};else,ZW=medium.WaveImpedance;end,Ts=(ZW^-1)*Ts;Ts=Ts*Z0;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "ABB"
                        Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(boundary.Dispersive),beta=boundary.Param(freqIndex);else,beta=boundary.Param;end,Ts=Ts*beta/k0;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "ABZ"
                       if(boundary.Tensor)
                            if(boundary.Dispersive),Z=boundary.Param{freqIndex};
                            else,Z=boundary.Param;
                            end
                            Ts=CalculateTensorIntegral(TModel,Z^-1,element,facets(kk),b,c,d,kk);Ts=Ts*Z0;
                        else
                            if(boundary.Dispersive),Z=boundary.Param(freqIndex);
                            else,Z=boundary.Param;
                            end
                            Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Ts=Ts*Z^-1;Ts=Ts*Z0;
                        end
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "GRA"
                        Tg = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                        if(boundary.Dispersive),cond=boundary.Param(freqIndex);else,cond=boundary.Param;end,Tg=Tg.*cond;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Tg(ii,jj)=Tg(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "IBC"
                        if(boundary.Dispersive),cond=boundary.Param(freqIndex);else,cond=boundary.Param;end,Tbc=Tbc.*cond;
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6,for jj=1:6,Tbc(ii,jj)=Tbc(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);end,end
                        end
                    case "POR"
                        if(boundary.Dispersive)
                            if(boundary.PortParamType==0)
                                Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Ts=Ts*boundary.Param(freqIndex)/k0;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==1),Ts=CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Zp=boundary.Param(freqIndex);
                                Ts=(Z0/Zp)*Ts;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==2)
                                Zp=boundary.Param{freqIndex};
                                Ts=CalculateTensorIntegral(TModel,Zp^-1,element,facets(kk),b,c,d,kk);
                                Ts=Z0*Ts;Tp=2*1i*Ts;
                            end
                        else
                            if(boundary.PortParamType==0)
                                Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);
                                Ts=boundary.Param*Ts/k0;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==1),Ts = CalculateIntegral(TModel,element,facets(kk),b,c,d,kk);Zp=boundary.Param; 
                                Ts=(Z0/Zp)*Ts;Tp=2*1i*Ts;
                            elseif(boundary.PortParamType==2),Zp=boundary.Param;Ts=CalculateTensorIntegral(TModel,Zp^-1,element,facets(kk),b,c,d,kk);
                                Ts=Z0*Ts;Tp=2*1i*Ts;
                            end
                        end
                        switch TModel.Assembled.E_Scaling
                            case 1,for ii=1:6,for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*edgeLengths(ii)*edgeLengths(jj);Ts(ii,jj)=Ts(ii,jj)*edgeLengths(ii)*edgeLengths(jj);end,end
                            case 2,for ii=1:6
                                        for jj=1:6,Tp(ii,jj)=Tp(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);
                                                   Ts(ii,jj)=Ts(ii,jj)*TModel.Assembled.ScalingVec_E(edges(ii).Index)*TModel.Assembled.ScalingVec_E(edges(jj).Index);
                                        end
                                    end
                        end
                end
            end
        end
        %==================== Global Matrix Assembly ======================
        for ii=1:6,di=abs(edges(ii).UknownIndex);dsi=sign(edges(ii).UknownIndex);si=edgeSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);
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
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;AA(counterA)=si*dsi*sj*dsj*Am(ii,jj);
                                   counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TC(counterA)=si*dsi*sj*dsj*Tcm(ii,jj);
                end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:6,dj=abs(edges(jj).UknownIndex);dsj=sign(edges(jj).UknownIndex);sj=edgeSigns(jj);
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;FF(counterA)=si*dsi*sj*dsj*Fm(ii,jj);end
            end
        end
        for ii=1:4,di=abs(facets(ii).UknownIndex);dsi=sign(facets(ii).UknownIndex);si=facetSigns(ii);
            for jj=1:4,dj=abs(facets(jj).UknownIndex);dsj=sign(facets(jj).UknownIndex);sj=facetSigns(jj);
                if(di~=0 && dj~=0),counterA=counterA+1;IA(counterA)=di;JA(counterA)=dj;TB(counterA)=si*sj*dsi*dsj*Tb(ii,jj);end
            end
        end
     end
end
%============ Surface Integral Calculations ===============================
function [n] = findNormaVector(TModel,element,facetIndex),v1=TModel.Vertices(element.Vertices(1));v2=TModel.Vertices(element.Vertices(2));v3=TModel.Vertices(element.Vertices(3));v4=TModel.Vertices(element.Vertices(4));
    switch facetIndex
        case 1,vec1=[v2.X-v1.X;v2.Y-v1.Y;v2.Z-v1.Z];vec2=[v3.X-v1.X;v3.Y-v1.Y;v3.Z-v1.Z];vec3=[v1.X-v4.X;v1.Y-v4.Y;v1.Z-v4.Z];
        case 2,vec1=[v3.X-v2.X;v3.Y-v2.Y;v3.Z-v2.Z];vec2=[v4.X-v2.X;v4.Y-v2.Y;v4.Z-v2.Z];vec3=[v2.X-v1.X;v2.Y-v1.Y;v2.Z-v1.Z];
        case 3,vec1=[v3.X-v1.X;v3.Y-v1.Y;v3.Z-v1.Z];vec2=[v4.X-v1.X;v4.Y-v1.Y;v4.Z-v1.Z];vec3=[v1.X-v2.X;v1.Y-v2.Y;v1.Z-v2.Z];
        case 4,vec1=[v2.X-v1.X;v2.Y-v1.Y;v2.Z-v1.Z];vec2=[v4.X-v1.X;v4.Y-v1.Y;v4.Z-v1.Z];vec3=[v1.X-v3.X;v1.Y-v3.Y;v1.Z-v3.Z];
    end,n=cross(vec1,vec2);n=n/norm(n);if(n'*vec3>0),n=-n;end
end
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
