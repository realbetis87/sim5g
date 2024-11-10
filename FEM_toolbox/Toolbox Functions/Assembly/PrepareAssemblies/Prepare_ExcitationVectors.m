%--------------------------------------------------------------------------
%{

%}
%--------------------------------------------------------------------------
function [TModel] = Prepare_ExcitationVectors(TModel)
    if(isempty(TModel.Solution)),TModel.Solution=Solution();end
    TModel.Solution.Type="Excitation";
    UMax=max([TModel.Facets.UknownIndex]);KMax=max([TModel.Facets.KnownIndex]);
    if(TModel.Frequency.NF==1)
        if(any([TModel.Boundaries.Type]=="DIR")),TModel.Solution.KnownExcitation=zeros(KMax,1);end
        if(any([TModel.Boundaries.Type]=="POR")),TModel.Solution.UknownExcitation=zeros(UMax,1);end
        TModel.Solution.ExcitationVector=zeros(UMax,1);TModel.Solution.SolutionVector=zeros(UMax,1);
    else
        if(any([TModel.Boundaries.Type]=="DIR")),TModel.Solution.KnownExcitation=zeros(KMax,TModel.Frequency.NF);end
        if(any([TModel.Boundaries.Type]=="POR")),TModel.Solution.UknownExcitation=zeros(UMax,TModel.Frequency.NF);end
        TModel.Solution.ExcitationVector=zeros(UMax,TModel.Frequency.NF);TModel.Solution.SolutionVector=zeros(UMax,TModel.Frequency.NF);
    end
    for ie = 1 : numel(TModel.Boundary_Excitations),excitation=TModel.Boundary_Excitations(ie);
        switch excitation.Type
            case "DIR"
                if(TModel.Frequency.NF==1)
                    switch excitation.Assembly.Equation
                        case 1
                                for ii=1:numel(excitation.Edges),edge=TModel.Edges(excitation.Edges(ii));
                                    if(edge.KnownIndex~=0 && edge.IndexE~=0)
                                        TModel.Solution.KnownExcitation(edge.KnownIndex)=excitation.Vector(edge.IndexE);
                                    end
                                end
                                for ii=1:numel(excitation.Facets),facet=TModel.Facets(excitation.Facets(ii));
                                    if(facet.KnownIndex~=0)
                                      %  TModel.Solution.KnownExcitation(facet.KnownIndex)=ReturnFluxFieldOutOfField_TFLF(TModel,excitation,facet,excitation.Vector);
                                    end
                                end
                        case 2
                                for ii=1:numel(excitation.Edges),edge=TModel.Edges(excitation.Edges(ii));
                                    if(edge.KnownIndex~=0 && edge.IndexE~=0)
                                        TModel.Solution.KnownExcitation(edge.KnownIndex)=excitation.Vector(edge.IndexE);
                                    end
                                end
                                for ii=1:numel(excitation.Facets),facet=TModel.Facets(excitation.Facets(ii));
                                    if(facet.KnownIndex~=0)
                                        TModel.Solution.KnownExcitation(facet.KnownIndex)=ReturnFluxFieldOutOfField_TFTF(TModel,excitation,facet,excitation.Vector);
                                    end
                                end
                        case 3
                                for ii=1:numel(excitation.Edges),edge=TModel.Edges(excitation.Edges(ii));
                                    if(edge.KnownIndex~=0 && edge.IndexE~=0)
                                        TModel.Solution.KnownExcitation(edge.KnownIndex)=excitation.Vector(edge.IndexE);
                                    end
                                end
                                for ii=1:numel(excitation.Facets),facet=TModel.Facets(excitation.Facets(ii));
                                    if(facet.KnownIndex~=0)
                                        TModel.Solution.KnownExcitation(facet.KnownIndex)=ReturnFluxFieldOutOfField_VWE(TModel,excitation,facet,excitation.Vector);
                                    end
                                end
                    end
                else
                    for ff= 1 : TModel.Frequency.NF
                        switch excitation.Assembly.Equation
                             case 1
                                for ii=1:numel(excitation.Edges),edge=TModel.Edges(excitation.Edges(ii));exc=excitation.Vector{ff};
                                    if(edge.KnownIndex~=0 && edge.IndexE~=0)
                                        TModel.Solution.KnownExcitation(edge.KnownIndex,ff)=exc(edge.IndexE);
                                    end
                                end
                                for ii=1:numel(excitation.Facets),facet=TModel.Facets(excitation.Facets(ii));
                                    if(facet.KnownIndex~=0)
                                        TModel.Solution.KnownExcitation(facet.KnownIndex,ff)=ReturnFluxFieldOutOfField_TFLF(TModel,excitation,facet,ff);
                                    end
                                end
                            case 2
                                for ii=1:numel(excitation.Edges),edge=TModel.Edges(excitation.Edges(ii));
                                    if(edge.KnownIndex~=0 && edge.IndexE~=0)
                                        TModel.Solution.KnownExcitation(edge.KnownIndex,ff)=excitation.Vector(edge.IndexE,ff);
                                    end
                                end
                                for ii=1:numel(excitation.Facets),facet=TModel.Facets(excitation.Facets(ii));
                                    if(facet.KnownIndex~=0)
                                        TModel.Solution.KnownExcitation(facet.KnownIndex,ff)=ReturnFluxFieldOutOfField_TFTF(TModel,excitation,facet,ff);
                                    end
                                end
                            case 3
                                for ii=1:numel(excitation.Edges),edge=TModel.Edges(excitation.Edges(ii));
                                    if(edge.KnownIndex~=0 && edge.IndexE~=0)
                                        TModel.Solution.KnownExcitation(edge.KnownIndex,ff)=excitation.Vector(edge.IndexE,ff);
                                    end
                                end
                                for ii=1:numel(excitation.Facets),facet=TModel.Facets(excitation.Facets(ii));
                                    if(facet.KnownIndex~=0)
                                        TModel.Solution.KnownExcitation(facet.KnownIndex,ff)=ReturnFluxFieldOutOfField_VWE(TModel,excitation,facet,ff);
                                    end
                                end
                         end
                    end
                end
            case "POR"
                if(TModel.Frequency.NF==1)
                    for ii=1:numel(excitation.Edges),edge=TModel.Edges(excitation.Edges(ii));
                        if(edge.IndexE~=0),TModel.Solution.UknownExcitation(edge.UknownIndex)=excitation.Vector(edge.IndexE);end
                    end
                else
                    for ff=1:TModel.Frequency.NF,exc=excitation.Vector{ff};
                            for ii=1:numel(excitation.Edges),edge=TModel.Edges(excitation.Edges(ii));
                                if(edge.IndexE~=0),TModel.Solution.UknownExcitation(edge.UknownIndex,ff)=exc(edge.IndexE);end
                            end
                    end
                end
        end
    end
    if (any([TModel.Boundaries.PortType]==1))
        for ib = 1 : numel(TModel.Boundaries),boundary=TModel.Boundaries(ib);
            if(boundary.Type=="POR" && boundary.PortType==1)
                for ie = 1 : numel (boundary.Edges),edge=TModel.Edges(boundary.Edges(ie));
                    v1=TModel.Vertices(edge.Vertices(1));v2=TModel.Vertices(edge.Vertices(2));
                    edgeVector=[v2.X - v1.X; v2.Y -v1.Y;v2.Z-v1.Z];
                    if TModel.Frequency.NF==1
                        if(edge.UknownIndex~=0),TModel.Solution.UknownExcitation(edge.UknownIndex)=edgeVector'*boundary.PlaneWave;end
                    else
                        for ff=1:TModel.Frequency.NF,if(edge.UknownIndex~=0),TModel.Solution.UknownExcitation(edge.UknownIndex,ff)=edgeVector'*boundary.PlaneWave;end,end
                    end
                end
            elseif(boundary.Type=="DIR" && boundary.PortType==1)
                 for ie = 1 : numel (boundary.Edges),edge=TModel.Edges(boundary.Edges(ie));
                    v1=TModel.Vertices(edge.Vertices(1));v2=TModel.Vertices(edge.Vertices(2));
                    edgeVector=[v2.X - v1.X; v2.Y -v1.Y;v2.Z-v1.Z];
                    if TModel.Frequency.NF==1
                        if(edge.KnownIndex~=0),TModel.Solution.KnownExcitation(edge.KnownIndex)=edgeVector'*boundary.PlaneWave;end
                    else
                        for ff=1:TModel.Frequency.NF,if(edge.KnownIndex~=0),TModel.Solution.KnownExcitation(edge.KnownIndex,ff)=edgeVector'*boundary.PlaneWave;end,end
                    end
                end
            end
        end
    end
end
%==========================================================================
function [Field] = ReturnFluxFieldOutOfField_TFLF(varargin)
    if(nargin==4),TModel=varargin{1};excitation=varargin{2};facet=varargin{4};Vector=excitation.Vector;
    elseif(nargin==5),TModel=varargin{1};excitation=varargin{2};facet=varargin{4};freqIndex=varargin{5};Vector=excitation.Vector{freqIndex};
    end
    boundary=TModel.Boundaries(excitation.BoundaryIndices);medium=facet.Medium2D;Field=0;
    if(medium.IsDispersive),if(medium.Type=="Iso"),mu=medium.Mu(freqIndex);else,mu=medium.Mu{freqIndex};end
    else,mu=medium.Mu;
    end
    vertices=[TModel.Vertices(facet.Vertices)];x=[Vertices.X];y=[vertices.Y];z=[vertices.Z];
    xx=element.Barycenter(1);yy=element.Barycenter(2);zz=element.Barycenter(3);
    edges=[TModel.Edges(facet.Edges)];edgeSigns=[facet.EdgeSigns];
    switch abs(boundary.Axis)
        case 1
             De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');
             bb(1)=(z(2)-z(3))/De;  cc(1)=(y(3)-y(2))/De;  aa(1)=(y(2)*z(3)-y(3)*z(2))/De;
             bb(2)=(z(3)-z(1))/De;  cc(2)=(y(1)-y(3))/De;  aa(2)=(y(3)*z(1)-y(1)*z(3))/De;
             bb(3)=(z(1)-z(2))/De;  cc(3)=(y(2)-y(1))/De;  aa(3)=(y(1)*z(2)-y(2)*z(1))/De;
             %----------------------------------------------------------
             zeta(1)=aa(1)+bb(1)*yy+cc(1)*zz;      wwy(1)=ll(1)*(zeta(1)*bb(2)-zeta(2)*bb(1));         wwz(1)=ll(1)*(zeta(1)*cc(2)-zeta(2)*cc(1));   
             zeta(2)=aa(2)+bb(2)*yy+cc(2)*zz;      wwy(2)=ll(2)*(zeta(2)*bb(3)-zeta(3)*bb(2));         wwz(2)=ll(2)*(zeta(2)*cc(3)-zeta(3)*cc(2));
             zeta(3)=aa(3)+bb(3)*yy+cc(3)*zz;      wwy(3)=ll(3)*(zeta(3)*bb(1)-zeta(1)*bb(3));         wwz(3)=ll(3)*(zeta(3)*cc(1)-zeta(1)*cc(3));
             %--------------------------------------------------------------
             switch medium.Type
                 case "Iso"
                     if(vertices(kk).IndexH~=0),Field=Field+m0*mu*zeta(kk)*Vector(vertices(kk).IndexH);end
                 case "Anis"
                     if(vertices(kk).IndexH~=0),Field=Field+m0*mu(1,1)*zeta(kk)*Vector(vertices(kk).IndexH);end
                     if(edges(kk).IndexH~=0),Field=Field+m0*mu(1,2)*wwy(kk)*edgeSigns(kk)*Vector(edges(kk).IndexH);
                                             Field=Field+m0*mu(1,3)*wwz(kk)*edgeSigns(kk)*Vector(edges(kk).IndexH);
                     end
             end
             if(sign(boundary.Axis)<0),Field=-Field*facet.Surface;else,Field=Field*facet.Surface;end
        case 2
             De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');
             bb(1)=(z(2)-z(3))/De;  cc(1)=(x(3)-x(2))/De;  aa(1)=(x(2)*z(3)-x(3)*z(2))/De;
             bb(2)=(z(3)-z(1))/De;  cc(2)=(x(1)-x(3))/De;  aa(2)=(x(3)*z(1)-x(1)*z(3))/De;
             bb(3)=(z(1)-z(2))/De;  cc(3)=(x(2)-x(1))/De;  aa(3)=(x(1)*z(2)-x(2)*z(1))/De;
             %----------------------------------------------------------
             zeta(1)=aa(1)+bb(1)*xx+cc(1)*zz;      wwx(1)=ll(1)*(zeta(1)*bb(2)-zeta(2)*bb(1));         wwz(1)=ll(1)*(zeta(1)*cc(2)-zeta(2)*cc(1));   
             zeta(2)=aa(2)+bb(2)*xx+cc(2)*zz;      wwx(2)=ll(2)*(zeta(2)*bb(3)-zeta(3)*bb(2));         wwz(2)=ll(2)*(zeta(2)*cc(3)-zeta(3)*cc(2));
             zeta(3)=aa(3)+bb(3)*xx+cc(3)*zz;      wwx(3)=ll(3)*(zeta(3)*bb(1)-zeta(1)*bb(3));         wwz(3)=ll(3)*(zeta(3)*cc(1)-zeta(1)*cc(3));
             %--------------------------------------------------------------
             switch medium.Type
                 case "Iso"
                     if(vertices(kk).IndexH~=0),Field=Field+m0*mu*zeta(kk)*Vector(vertices(kk).IndexH);end
                 case "Anis"
                     if(vertices(kk).IndexH~=0),Field=Field+m0*mu(2,2)*zeta(kk)*Vector(vertices(kk).IndexH);end
                     if(edges(kk).IndexH~=0),Field=Field+m0*mu(2,1)*wwx(kk)*edgeSigns(kk)*Vector(edges(kk).IndexH);
                                             Field=Field+m0*mu(2,3)*wwz(kk)*edgeSigns(kk)*Vector(edges(kk).IndexH);
                     end
             end
             if(sign(boundary.Axis)<0),Field=-Field*facet.Surface;else,Field=Field*facet.Surface;end
        case 3
             De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');
             bb(1)=(y(2)-y(3))/De;  cc(1)=(x(3)-x(2))/De;  aa(1)=(x(2)*y(3)-x(3)*y(2))/De;
             bb(2)=(y(3)-y(1))/De;  cc(2)=(x(1)-x(3))/De;  aa(2)=(x(3)*y(1)-x(1)*y(3))/De;
             bb(3)=(y(1)-y(2))/De;  cc(3)=(x(2)-x(1))/De;  aa(3)=(x(1)*y(2)-x(2)*y(1))/De;
             %----------------------------------------------------------
             zeta(1)=aa(1)+bb(1)*xx+cc(1)*yy;      wwx(1)=ll(1)*(zeta(1)*bb(2)-zeta(2)*bb(1));         wwy(1)=ll(1)*(zeta(1)*cc(2)-zeta(2)*cc(1));   
             zeta(2)=aa(2)+bb(2)*xx+cc(2)*yy;      wwx(2)=ll(2)*(zeta(2)*bb(3)-zeta(3)*bb(2));         wwy(2)=ll(2)*(zeta(2)*cc(3)-zeta(3)*cc(2));
             zeta(3)=aa(3)+bb(3)*xx+cc(3)*yy;      wwx(3)=ll(3)*(zeta(3)*bb(1)-zeta(1)*bb(3));         wwy(3)=ll(3)*(zeta(3)*cc(1)-zeta(1)*cc(3));
             %--------------------------------------------------------------
             switch medium.Type
                 case "Iso"
                     if(vertices(kk).IndexH~=0),Field=Field+m0*mu*zeta(kk)*Vector(vertices(kk).IndexH);end
                 case "Anis"
                     if(vertices(kk).IndexH~=0),Field=Field+m0*mu(3,3)*zeta(kk)*Vector(vertices(kk).IndexH);end
                     if(edges(kk).IndexH~=0),Field=Field+m0*mu(3,1)*wwx(kk)*edgeSigns(kk)*Vector(edges(kk).IndexH);
                                             Field=Field+m0*mu(3,2)*wwy(kk)*edgeSigns(kk)*Vector(edges(kk).IndexH);
                     end
             end
             if(sign(boundary.Axis)<0),Field=-Field*facet.Surface;else,Field=Field*facet.Surface;end
     
    end
end
