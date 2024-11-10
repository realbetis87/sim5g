function [Field,xx,yy,zz] = PlotFieldOnSlice(XU,XK,FieldName,FieldQ,Position,Points),global Edges;global Facets;global NumberOfFacets;global NumberOfEdges;global ComputationalDomain;Constants;
%--------------------------------------------------------------------------
%{
       1. XU        : Uknown Field to be plotted (1 x Total Number Of  Uknowns)
       2. XK        : Known Field to be plotted  (1 x Total Number Of  Knowns)
       2. Position  : Plane Coordinates - for (-Ax) plane
                                Position=[-Ax;Nan;Nan;]
       3.Points     :[NumberOfPoints 1 Axis, NumberOfPOints 2 Axis]
%}
%--------------------------------------------------------------------------
    FacetVector=zeros(NumberOfFacets,1);EdgeVector=zeros(NumberOfEdges,1);
    for ii=1:NumberOfFacets,facet=Facets(ii);if(facet.KnownIndex~=0),FacetVector(ii)=XK(facet.KnownIndex);elseif(facet.UknownIndex~=0),FacetVector(ii)=XU(facet.UknownIndex);end,end,FacetVector=-(1i/c0)*FacetVector;
    for ii=1:NumberOfEdges,edge=Edges(ii);if(edge.KnownIndex~=0),EdgeVector(ii)=XK(edge.KnownIndex);elseif(edge.UknownIndex~=0),EdgeVector(ii)=XU(edge.UknownIndex);end,end
    Plane=find(~isnan(Position));N1=Points(1);N2=Points(2);empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    Xmin=min(ComputationalDomain.XCorners);Xmax=max(ComputationalDomain.XCorners);Ymin=min(ComputationalDomain.YCorners);Ymax=max(ComputationalDomain.YCorners);Zmin=min(ComputationalDomain.ZCorners);Zmax=max(ComputationalDomain.ZCorners);
    Wx=ComputationalDomain.XDimension;Wy=ComputationalDomain.YDimension;Wz=ComputationalDomain.ZDimension;
    Xmin=Xmin+Wx/100;Xmax=Xmax-Wx/100;Ymin=Ymin+Wy/100;Ymax=Ymax-Wy/100;Zmin=Zmin+Wz/100;Zmax=Zmax-Wz/100;
    switch Plane
        case 1,xx=Position(1);yy=linspace(Ymin,Ymax,N1);zz=linspace(Zmin,Zmax,N2);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotFieldXPlane(FieldName,FieldQ,EdgeVector,FacetVector,xx,yy,zz,N1,N2);
        case 2,xx=linspace(Xmin,Xmax,N1);yy=Position(2);zz=linspace(Zmin,Zmax,N2);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotFieldYPlane(FieldName,FieldQ,EdgeVector,FacetVector,xx,yy,zz,N1,N2);
        case 3,xx=linspace(Xmin,Xmax,N1);yy=linspace(Ymin,Ymax,N2);zz=Position(3);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotFieldZPlane(FieldName,FieldQ,EdgeVector,FacetVector,xx,yy,zz,N1,N2);
    end
    switch FieldName
        case "Ex",Field=Ex;
        case "Ey",Field=Ey;
        case "Ez",Field=Ez;
        case "Bx",Field=Bx;
        case "By",Field=By;
        case "Bz",Field=Bz;
    end
end
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotFieldXPlane(FieldName,FieldQ,EdgeVector,FacetVector,xx,yy,zz,N1,N2),global Elements;global Edges;global Facets;global error;empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    Barycenter=[Elements.Barycenter];disp("X Plane Plot")
    for ii=1:N1
        for jj=1:N2,distance=((xx-Barycenter(1,:)).^2+(yy(ii)-Barycenter(2,:)).^2 + (zz(jj)-Barycenter(3,:)).^2);[~,Is]=sort(distance);ie=0;aux=0;
            while(ie==0),aux=aux+1;element=Elements(Is(aux));a=element.As;b=element.Bs;c=element.Cs;d=element.Ds;
                z(1)=a(1)+b(1)*xx+c(1)*yy(ii)+d(1)*zz(jj);z(2)=a(2)+b(2)*xx+c(2)*yy(ii)+d(2)*zz(jj);z(3)=a(3)+b(3)*xx+c(3)*yy(ii)+d(3)*zz(jj);z(4)=a(4)+b(4)*xx+c(4)*yy(ii)+d(4)*zz(jj);
                if (z(1)>0 || abs(z(1))<=error) && (z(2)>0 || abs(z(2))<=error) && (z(3)>0 || abs(z(3))<=error)&& (z(4)>0 || abs(z(4))<=error),ie=Is(aux);end
            end,edges=[element.Edges];facets=[element.Facets];ll=[Edges(edges).Length];ss=[Facets(facets).Surface];
            wfx(1)=2*z(3)*(c(1)*d(2)-c(2)*d(1))+2*z(1)*(c(2)*d(3)-c(3)*d(2))+2*z(2)*(c(3)*d(1)-c(1)*d(3));wfy(1)=2*z(3)*(d(1)*b(2)-d(2)*b(1))+2*z(1)*(d(2)*b(3)-d(3)*b(2))+2*z(2)*(d(3)*b(1)-d(1)*b(3));wfz(1)=2*z(3)*(b(1)*c(2)-b(2)*c(1))+2*z(1)*(b(2)*c(3)-b(3)*c(2))+2*z(2)*(b(3)*c(1)-b(1)*c(3));
            wfx(2)=2*z(3)*(c(2)*d(4)-c(4)*d(2))+2*z(2)*(c(4)*d(3)-c(3)*d(4))+2*z(4)*(c(3)*d(2)-c(2)*d(3));wfy(2)=2*z(3)*(d(2)*b(4)-d(4)*b(2))+2*z(2)*(d(4)*b(3)-d(3)*b(4))+2*z(4)*(d(3)*b(2)-d(2)*b(3));wfz(2)=2*z(3)*(b(2)*c(4)-b(4)*c(2))+2*z(2)*(b(4)*c(3)-b(3)*c(4))+2*z(4)*(b(3)*c(2)-b(2)*c(3));
            wfx(3)=2*z(1)*(c(3)*d(4)-c(4)*d(3))+2*z(2)*(c(4)*d(1)-c(1)*d(4))+2*z(4)*(c(1)*d(3)-c(3)*d(1));wfy(3)=2*z(1)*(d(3)*b(4)-d(4)*b(3))+2*z(2)*(d(4)*b(1)-d(1)*b(4))+2*z(4)*(d(1)*b(3)-d(3)*b(1));wfz(3)=2*z(1)*(b(3)*c(4)-b(4)*c(3))+2*z(2)*(b(4)*c(1)-b(1)*c(4))+2*z(4)*(b(1)*c(3)-b(3)*c(1));
            wfx(4)=2*z(1)*(c(4)*d(2)-c(2)*d(4))+2*z(4)*(c(2)*d(1)-c(1)*d(2))+2*z(2)*(c(1)*d(4)-c(4)*d(1));wfy(4)=2*z(1)*(d(4)*b(2)-d(2)*b(4))+2*z(4)*(d(2)*b(1)-d(1)*b(2))+2*z(2)*(d(1)*b(4)-d(4)*b(1));wfz(4)=2*z(1)*(b(4)*c(2)-b(2)*c(4))+2*z(4)*(b(2)*c(1)-b(1)*c(2))+2*z(2)*(b(1)*c(4)-b(4)*c(1));
            wx(1)=z(1)*b(2)-z(2)*b(1); wx(2)=z(1)*b(3)-z(3)*b(1);wx(3)=z(1)*b(4)-z(4)*b(1);wx(4)=z(2)*b(3)-z(3)*b(2);wx(5)=z(2)*b(4)-z(4)*b(2);wx(6)=z(3)*b(4)-z(4)*b(3);wy(1)=z(1)*c(2)-z(2)*c(1); wy(2)=z(1)*c(3)-z(3)*c(1);wy(3)=z(1)*c(4)-z(4)*c(1);wy(4)=z(2)*c(3)-z(3)*c(2);wy(5)=z(2)*c(4)-z(4)*c(2);wy(6)=z(3)*c(4)-z(4)*c(3);
            wz(1)=z(1)*d(2)-z(2)*d(1); wz(2)=z(1)*d(3)-z(3)*d(1);wz(3)=z(1)*d(4)-z(4)*d(1);wz(4)=z(2)*d(3)-z(3)*d(2);wz(5)=z(2)*d(4)-z(4)*d(2);wz(6)=z(3)*d(4)-z(4)*d(3);
            wx=wx.*ll;wy=wy.*ll;wz=wz.*ll;wfx=wfx.*ss;wfy=wfy.*ss;wfz=wfz.*ss;
            for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*sgn*FacetVector(facets(kk));By(ii,jj)=By(ii,jj)+wfy(kk)*sgn*FacetVector(facets(kk));Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*sgn*FacetVector(facets(kk));end
            for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*sgn*EdgeVector(edges(kk));Ey(ii,jj)=Ey(ii,jj)+wy(kk)*sgn*EdgeVector(edges(kk));Ez(ii,jj)=Ez(ii,jj)+wz(kk)*sgn*EdgeVector(edges(kk));end
        end
     end
    [cmap]=buildcmap('cbkry');figure;
    switch FieldName
        case "Ex"
            switch FieldQ
                case "r",pcolor(real(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field absolute');colorbar;axis tight;axis equal;
            end
        case "Ey"
            switch FieldQ
                case "r",pcolor(real(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field absolute');colorbar;axis tight;axis equal;
            end
        case "Ez"
            switch FieldQ
                case "r",pcolor(real(Ez));shading interp;hold on;colormap(cmap);title('z Component Electric Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Ez));shading interp;hold on;colormap(cmap);title('z Component Electric Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Ez));shading interp;hold on;colormap(cmap);title('z Component Electric Field absolute');colorbar;axis tight;axis equal;
            end
        case "Bx"
            switch FieldQ
                case "r",pcolor(real(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field absolute');colorbar;axis tight;axis equal;
            end
        case "By"
            switch FieldQ
                case "r",pcolor(real(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field absolute');colorbar;axis tight;axis equal;
            end
        case "Bz"
            switch FieldQ
                case "r",pcolor(real(Bz));shading interp;hold on;colormap(cmap);title('z Component Magnetic Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Bz));shading interp;hold on;colormap(cmap);title('z Component Magnetic Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Bz));shading interp;hold on;colormap(cmap);title('z Component Magnetic Field absolute');colorbar;axis tight;axis equal;
            end
    end
end
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotFieldYPlane(FieldName,FieldQ,EdgeVector,FacetVector,xx,yy,zz,N1,N2),global Elements;global Edges;global Facets;global error;empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    Barycenter=[Elements.Barycenter];disp("Z Plane Plot")
    for ii=1:N1
        for jj=1:N2,distance=((xx(ii)-Barycenter(1,:)).^2+(yy-Barycenter(2,:)).^2 + (zz(jj)-Barycenter(3,:)).^2);[~,Is]=sort(distance);ie=0;aux=0;
            while(ie==0),aux=aux+1;element=Elements(Is(aux));a=element.As;b=element.Bs;c=element.Cs;d=element.Ds;
                z(1)=a(1)+b(1)*xx(ii)+c(1)*yy+d(1)*zz(jj);z(2)=a(2)+b(2)*xx(ii)+c(2)*yy+d(2)*zz(jj);z(3)=a(3)+b(3)*xx(ii)+c(3)*yy+d(3)*zz(jj);z(4)=a(4)+b(4)*xx(ii)+c(4)*yy+d(4)*zz(jj);
                if (z(1)>0 || abs(z(1))<=error) && (z(2)>0 || abs(z(2))<=error) && (z(3)>0 || abs(z(3))<=error)&& (z(4)>0 || abs(z(4))<=error),ie=Is(aux);end
            end,edges=[element.Edges];facets=[element.Facets];ll=[Edges(edges).Length];ss=[Facets(facets).Surface];
            wfx(1)=2*z(3)*(c(1)*d(2)-c(2)*d(1))+2*z(1)*(c(2)*d(3)-c(3)*d(2))+2*z(2)*(c(3)*d(1)-c(1)*d(3));wfy(1)=2*z(3)*(d(1)*b(2)-d(2)*b(1))+2*z(1)*(d(2)*b(3)-d(3)*b(2))+2*z(2)*(d(3)*b(1)-d(1)*b(3));wfz(1)=2*z(3)*(b(1)*c(2)-b(2)*c(1))+2*z(1)*(b(2)*c(3)-b(3)*c(2))+2*z(2)*(b(3)*c(1)-b(1)*c(3));
            wfx(2)=2*z(3)*(c(2)*d(4)-c(4)*d(2))+2*z(2)*(c(4)*d(3)-c(3)*d(4))+2*z(4)*(c(3)*d(2)-c(2)*d(3));wfy(2)=2*z(3)*(d(2)*b(4)-d(4)*b(2))+2*z(2)*(d(4)*b(3)-d(3)*b(4))+2*z(4)*(d(3)*b(2)-d(2)*b(3));wfz(2)=2*z(3)*(b(2)*c(4)-b(4)*c(2))+2*z(2)*(b(4)*c(3)-b(3)*c(4))+2*z(4)*(b(3)*c(2)-b(2)*c(3));
            wfx(3)=2*z(1)*(c(3)*d(4)-c(4)*d(3))+2*z(2)*(c(4)*d(1)-c(1)*d(4))+2*z(4)*(c(1)*d(3)-c(3)*d(1));wfy(3)=2*z(1)*(d(3)*b(4)-d(4)*b(3))+2*z(2)*(d(4)*b(1)-d(1)*b(4))+2*z(4)*(d(1)*b(3)-d(3)*b(1));wfz(3)=2*z(1)*(b(3)*c(4)-b(4)*c(3))+2*z(2)*(b(4)*c(1)-b(1)*c(4))+2*z(4)*(b(1)*c(3)-b(3)*c(1));
            wfx(4)=2*z(1)*(c(4)*d(2)-c(2)*d(4))+2*z(4)*(c(2)*d(1)-c(1)*d(2))+2*z(2)*(c(1)*d(4)-c(4)*d(1));wfy(4)=2*z(1)*(d(4)*b(2)-d(2)*b(4))+2*z(4)*(d(2)*b(1)-d(1)*b(2))+2*z(2)*(d(1)*b(4)-d(4)*b(1));wfz(4)=2*z(1)*(b(4)*c(2)-b(2)*c(4))+2*z(4)*(b(2)*c(1)-b(1)*c(2))+2*z(2)*(b(1)*c(4)-b(4)*c(1));
            wx(1)=z(1)*b(2)-z(2)*b(1); wx(2)=z(1)*b(3)-z(3)*b(1);wx(3)=z(1)*b(4)-z(4)*b(1);wx(4)=z(2)*b(3)-z(3)*b(2);wx(5)=z(2)*b(4)-z(4)*b(2);wx(6)=z(3)*b(4)-z(4)*b(3);wy(1)=z(1)*c(2)-z(2)*c(1); wy(2)=z(1)*c(3)-z(3)*c(1);wy(3)=z(1)*c(4)-z(4)*c(1);wy(4)=z(2)*c(3)-z(3)*c(2);wy(5)=z(2)*c(4)-z(4)*c(2);wy(6)=z(3)*c(4)-z(4)*c(3);
            wz(1)=z(1)*d(2)-z(2)*d(1); wz(2)=z(1)*d(3)-z(3)*d(1);wz(3)=z(1)*d(4)-z(4)*d(1);wz(4)=z(2)*d(3)-z(3)*d(2);wz(5)=z(2)*d(4)-z(4)*d(2);wz(6)=z(3)*d(4)-z(4)*d(3);
            wx=wx.*ll;wy=wy.*ll;wz=wz.*ll;wfx=wfx.*ss;wfy=wfy.*ss;wfz=wfz.*ss;
            for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*sgn*FacetVector(facets(kk));By(ii,jj)=By(ii,jj)+wfy(kk)*sgn*FacetVector(facets(kk));Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*sgn*FacetVector(facets(kk));end
            for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*sgn*EdgeVector(edges(kk));Ey(ii,jj)=Ey(ii,jj)+wy(kk)*sgn*EdgeVector(edges(kk));Ez(ii,jj)=Ez(ii,jj)+wz(kk)*sgn*EdgeVector(edges(kk));
                    if(real(Ez(ii,jj))<0),Ex(ii,jj)=0;end
            end
        end
    end
   [cmap]=buildcmap('cbkry');figure;
    switch FieldName
        case "Ex"
            switch FieldQ
                case "r",pcolor(real(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field absolute');colorbar;axis tight;axis equal;
            end
        case "Ey"
            switch FieldQ
                case "r",pcolor(real(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field absolute');colorbar;axis tight;axis equal;
            end
        case "Ez"
            switch FieldQ
                case "r",pcolor(real(Ez));shading interp;hold on;colormap(cmap);title('z Component Electric Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Ez));shading interp;hold on;colormap(cmap);title('z Component Electric Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Ez));shading interp;hold on;colormap(cmap);title('z Component Electric Field absolute');colorbar;axis tight;axis equal;
            end
        case "Bx"
            switch FieldQ
                case "r",pcolor(real(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field absolute');colorbar;axis tight;axis equal;
            end
        case "By"
            switch FieldQ
                case "r",pcolor(real(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field absolute');colorbar;axis tight;axis equal;
            end
        case "Bz"
            switch FieldQ
                case "r",pcolor(real(Bz));shading interp;hold on;colormap(cmap);title('z Component Magnetic Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Bz));shading interp;hold on;colormap(cmap);title('z Component Magnetic Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Bz));shading interp;hold on;colormap(cmap);title('z Component Magnetic Field absolute');colorbar;axis tight;axis equal;
            end
    end
end
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotFieldZPlane(FieldName,FieldQ,EdgeVector,FacetVector,xx,yy,zz,N1,N2),global Elements;global Edges;global Facets;global error;empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    Barycenter=[Elements.Barycenter];disp("Z Plane Plot")
    for ii=1:N1
        for jj=1:N2,distance=((xx(ii)-Barycenter(1,:)).^2+(yy(jj)-Barycenter(2,:)).^2 + (zz-Barycenter(3,:)).^2);[~,Is]=sort(distance);ie=0;aux=0;
            while(ie==0),aux=aux+1;element=Elements(Is(aux));a=element.As;b=element.Bs;c=element.Cs;d=element.Ds;
                z(1)=a(1)+b(1)*xx(ii)+c(1)*yy(jj)+d(1)*zz;z(2)=a(2)+b(2)*xx(ii)+c(2)*yy(jj)+d(2)*zz;z(3)=a(3)+b(3)*xx(ii)+c(3)*yy(jj)+d(3)*zz;z(4)=a(4)+b(4)*xx(ii)+c(4)*yy(jj)+d(4)*zz;
                if (z(1)>0 || abs(z(1))<=error) && (z(2)>0 || abs(z(2))<=error) && (z(3)>0 || abs(z(3))<=error)&& (z(4)>0 || abs(z(4))<=error),ie=Is(aux);end
            end,edges=[element.Edges];facets=[element.Facets];ll=[Edges(edges).Length];ss=[Facets(facets).Surface];
            wfx(1)=2*z(3)*(c(1)*d(2)-c(2)*d(1))+2*z(1)*(c(2)*d(3)-c(3)*d(2))+2*z(2)*(c(3)*d(1)-c(1)*d(3));wfy(1)=2*z(3)*(d(1)*b(2)-d(2)*b(1))+2*z(1)*(d(2)*b(3)-d(3)*b(2))+2*z(2)*(d(3)*b(1)-d(1)*b(3));wfz(1)=2*z(3)*(b(1)*c(2)-b(2)*c(1))+2*z(1)*(b(2)*c(3)-b(3)*c(2))+2*z(2)*(b(3)*c(1)-b(1)*c(3));
            wfx(2)=2*z(3)*(c(2)*d(4)-c(4)*d(2))+2*z(2)*(c(4)*d(3)-c(3)*d(4))+2*z(4)*(c(3)*d(2)-c(2)*d(3));wfy(2)=2*z(3)*(d(2)*b(4)-d(4)*b(2))+2*z(2)*(d(4)*b(3)-d(3)*b(4))+2*z(4)*(d(3)*b(2)-d(2)*b(3));wfz(2)=2*z(3)*(b(2)*c(4)-b(4)*c(2))+2*z(2)*(b(4)*c(3)-b(3)*c(4))+2*z(4)*(b(3)*c(2)-b(2)*c(3));
            wfx(3)=2*z(1)*(c(3)*d(4)-c(4)*d(3))+2*z(2)*(c(4)*d(1)-c(1)*d(4))+2*z(4)*(c(1)*d(3)-c(3)*d(1));wfy(3)=2*z(1)*(d(3)*b(4)-d(4)*b(3))+2*z(2)*(d(4)*b(1)-d(1)*b(4))+2*z(4)*(d(1)*b(3)-d(3)*b(1));wfz(3)=2*z(1)*(b(3)*c(4)-b(4)*c(3))+2*z(2)*(b(4)*c(1)-b(1)*c(4))+2*z(4)*(b(1)*c(3)-b(3)*c(1));
            wfx(4)=2*z(1)*(c(4)*d(2)-c(2)*d(4))+2*z(4)*(c(2)*d(1)-c(1)*d(2))+2*z(2)*(c(1)*d(4)-c(4)*d(1));wfy(4)=2*z(1)*(d(4)*b(2)-d(2)*b(4))+2*z(4)*(d(2)*b(1)-d(1)*b(2))+2*z(2)*(d(1)*b(4)-d(4)*b(1));wfz(4)=2*z(1)*(b(4)*c(2)-b(2)*c(4))+2*z(4)*(b(2)*c(1)-b(1)*c(2))+2*z(2)*(b(1)*c(4)-b(4)*c(1));
            wx(1)=z(1)*b(2)-z(2)*b(1); wx(2)=z(1)*b(3)-z(3)*b(1);wx(3)=z(1)*b(4)-z(4)*b(1);wx(4)=z(2)*b(3)-z(3)*b(2);wx(5)=z(2)*b(4)-z(4)*b(2);wx(6)=z(3)*b(4)-z(4)*b(3);wy(1)=z(1)*c(2)-z(2)*c(1); wy(2)=z(1)*c(3)-z(3)*c(1);wy(3)=z(1)*c(4)-z(4)*c(1);wy(4)=z(2)*c(3)-z(3)*c(2);wy(5)=z(2)*c(4)-z(4)*c(2);wy(6)=z(3)*c(4)-z(4)*c(3);
            wz(1)=z(1)*d(2)-z(2)*d(1); wz(2)=z(1)*d(3)-z(3)*d(1);wz(3)=z(1)*d(4)-z(4)*d(1);wz(4)=z(2)*d(3)-z(3)*d(2);wz(5)=z(2)*d(4)-z(4)*d(2);wz(6)=z(3)*d(4)-z(4)*d(3);
            wx=wx.*ll;wy=wy.*ll;wz=wz.*ll;wfx=wfx.*ss;wfy=wfy.*ss;wfz=wfz.*ss;
            for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*sgn*FacetVector(facets(kk));By(ii,jj)=By(ii,jj)+wfy(kk)*sgn*FacetVector(facets(kk));Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*sgn*FacetVector(facets(kk));end
            for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*sgn*EdgeVector(edges(kk));Ey(ii,jj)=Ey(ii,jj)+wy(kk)*sgn*EdgeVector(edges(kk));Ez(ii,jj)=Ez(ii,jj)+wz(kk)*sgn*EdgeVector(edges(kk));end
        end
    end,Bx=Bx(:,35:end);disp("here");
    [cmap]=buildcmap('cbkry');figure;
    switch FieldName
        case "Ex"
            switch FieldQ
                case "r",pcolor(real(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field absolute');colorbar;axis tight;axis equal;
            end
        case "Ey"
            switch FieldQ
                case "r",pcolor(real(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field absolute');colorbar;axis tight;axis equal;
            end
        case "Ez"
            switch FieldQ
                case "r",pcolor(real(Ez));shading interp;hold on;colormap(cmap);title('z Component Electric Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Ez));shading interp;hold on;colormap(cmap);title('z Component Electric Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Ez));shading interp;hold on;colormap(cmap);title('z Component Electric Field absolute');colorbar;axis tight;axis equal;
            end
        case "Bx"
            switch FieldQ
                case "r",pcolor(real(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field absolute');colorbar;axis tight;axis equal;
            end
        case "By"
            switch FieldQ
                case "r",pcolor(real(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field absolute');colorbar;axis tight;axis equal;
            end
        case "Bz"
            switch FieldQ
                case "r",pcolor(real(Bz));shading interp;hold on;colormap(cmap);title('z Component Magnetic Field real');colorbar;axis tight;axis equal;
                case "i",pcolor(imag(Bz));shading interp;hold on;colormap(cmap);title('z Component Magnetic Field imaginary');colorbar;axis tight;axis equal;
                case "a",pcolor(abs(Bz));shading interp;hold on;colormap(cmap);title('z Component Magnetic Field absolute');colorbar;axis tight;axis equal;
            end
    end
end