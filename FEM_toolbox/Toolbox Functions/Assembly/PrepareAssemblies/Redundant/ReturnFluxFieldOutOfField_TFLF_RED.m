function [Field] = ReturnFluxFieldOutOfField_TFLF(varargin)
    if(nargin==4),TModel=varargin{1};excitation=varargin{2};Vector=varargin{3};facet=varargin{4};
    elseif(nargin==5),TModel=varargin{1};excitation=varargin{2};Vector=varargin{3};facet=varargin{4};freqIndex=varargin{5};
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