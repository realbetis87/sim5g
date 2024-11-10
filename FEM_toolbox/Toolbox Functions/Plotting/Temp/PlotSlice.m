function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz] = PlotSlice(XU,XK,Position,Points),global Edges;global Facets;global NumberOfFacets;global NumberOfEdges;global ComputationalDomain;Constants;
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
        case 1,xx=Position(1);yy=linspace(Ymin,Ymax,N1);zz=linspace(Zmin,Zmax,N2);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotXPlane(EdgeVector,FacetVector,xx,yy,zz,N1,N2);
        case 2,xx=linspace(Xmin,Xmax,N1);yy=Position(2);zz=linspace(Zmin,Zmax,N2);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotYPlane(EdgeVector,FacetVector,xx,yy,zz,N1,N2);
        case 3,xx=linspace(Xmin,Xmax,N1);yy=linspace(Ymin,Ymax,N2);zz=Position(3);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotZPlane(EdgeVector,FacetVector,xx,yy,zz,N1,N2);
    end
end

function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotXPlane(EdgeVector,FacetVector,xx,yy,zz,N1,N2),global Elements;global Edges;global Facets;global error;empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    Barycenter=[Elements.Barycenter];
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
    subplot(2,3,1);pcolor(real(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(real(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(real(Ez));shading interp;colormap(cmap);title('z Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(real(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(real(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(real(Bz));shading interp;colormap(cmap);title('z Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(2,3,1);pcolor(imag(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(imag(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(imag(Ez));shading interp;colormap(cmap);title('z Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(imag(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(imag(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(imag(Bz));shading interp;colormap(cmap);title('z Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
end
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotYPlane(EdgeVector,FacetVector,xx,yy,zz,N1,N2),global Elements;global Edges;global Facets;global error;empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    Barycenter=[Elements.Barycenter];
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
            for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*sgn*EdgeVector(edges(kk));Ey(ii,jj)=Ey(ii,jj)+wy(kk)*sgn*EdgeVector(edges(kk));Ez(ii,jj)=Ez(ii,jj)+wz(kk)*sgn*EdgeVector(edges(kk));end
        end
     end
    [cmap]=buildcmap('cbkry');figure;
    subplot(2,3,1);pcolor(real(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(real(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(real(Ez));shading interp;colormap(cmap);title('z Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(real(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(real(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(real(Bz));shading interp;colormap(cmap);title('z Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(2,3,1);pcolor(imag(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(imag(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(imag(Ez));shading interp;colormap(cmap);title('z Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(imag(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(imag(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(imag(Bz));shading interp;colormap(cmap);title('z Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
end
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotZPlane(EdgeVector,FacetVector,xx,yy,zz,N1,N2),global Elements;global Edges;global Facets;global error;empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    Barycenter=[Elements.Barycenter];
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
     end
    [cmap]=buildcmap('cbkry');figure;
    subplot(2,3,1);pcolor(real(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(real(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(real(Ez));shading interp;colormap(cmap);title('z Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(real(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(real(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(real(Bz));shading interp;colormap(cmap);title('z Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(2,3,1);pcolor(imag(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(imag(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(imag(Ez));shading interp;colormap(cmap);title('z Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(imag(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(imag(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(imag(Bz));shading interp;colormap(cmap);title('z Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
end

%{
    global Vertices;global Elements;global Facets;global Edges;global error;global ComputationalDomain;global NumberOfFacets;global NumberOfEdges;Constants;
    FacetVector=zeros(NumberOfFacets,1);EdgeVector=zeros(NumberOfEdges,1);
    for ii=1:NumberOfFacets,facet=Facets(ii);if(facet.KnownIndex~=0),FacetVector(ii)=XK(facet.KnownIndex);elseif(facet.UknownIndex~=0),FacetVector(ii)=XU(facet.UknownIndex);end,end%,FacetVector=-(1i/c0)*FacetVector;
    for ii=1:NumberOfEdges,edge=Edges(ii);if(edge.KnownIndex~=0),EdgeVector(ii)=XK(edge.KnownIndex);elseif(edge.UknownIndex~=0),EdgeVector(ii)=XU(edge.UknownIndex);end,end
    Plane=find(~isnan(Position));N1=Points(1);N2=Points(2);empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;Barycenter=[Elements.Barycenter];
    Xmin=min(ComputationalDomain.XCorners);Xmax=max(ComputationalDomain.XCorners);Ymin=min(ComputationalDomain.YCorners);Ymax=max(ComputationalDomain.YCorners);Zmin=min(ComputationalDomain.ZCorners);Zmax=max(ComputationalDomain.ZCorners);
    Wx=ComputationalDomain.XDimension;Wy=ComputationalDomain.YDimension;Wz=ComputationalDomain.ZDimension;
    Xmin=Xmin+Wx/5;Xmax=Xmax-Wx/100;Ymin=Ymin+Wy/100;Ymax=Ymax-Wy/100;Zmin=Zmin+Wz/100;Zmax=Zmax-Wz/100;
    switch Plane
        case 1,xx=Position(1);yy=linspace(Ymin,Ymax,N1);zz=linspace(Zmin,Zmax,N2);L1=yy;L2=zz;
        case 2,xx=linspace(Xmin,Xmax,N1);yy=Position(2);zz=linspace(Zmin,Zmax,N2);L1=xx;L2=zz;
        case 3,xx=linspace(Xmin,Xmax,N1);yy=linspace(Ymin,Ymax,N2);zz=Position(3);L1=yy;L2=xx;
    end
    for ii=1:N1
        for jj=1:N2
            switch Plane
                case 1,distance=((xx-Barycenter(1,:)).^2+(yy(ii)-Barycenter(2,:)).^2 + (zz(jj)-Barycenter(3,:)).^2);
                case 2,distance=((xx(ii)-Barycenter(1,:)).^2+(yy-Barycenter(2,:)).^2 + (zz(jj)-Barycenter(3,:)).^2);
                case 3,distance=((xx(ii)-Barycenter(1,:)).^2 +(yy(jj)-Barycenter(2,:)).^2+(zz-Barycenter(3,:)).^2);
            end,[~,Is]=sort(distance);ie=0;aux=0;
            while(ie==0),aux=aux+1;element=Elements(Is(aux));Nodes=element.Vertices;x=[Vertices(Nodes).X];y=[Vertices(Nodes).Y];zd=[Vertices(Nodes).Z];edges=[element.Edges];facets=[element.Facets];
                  De=det([1 x(1) y(1) zd(1);1 x(2) y(2) zd(2);1 x(3) y(3) zd(3);1 x(4) y(4) zd(4);]);ll=[Edges(edges).Length];ss=[Facets(facets).Surface];%ll=ones(6,1);
                  a(1)=det ([1 x(1) y(1) zd(1); 0 x(2) y(2) zd(2); 0 x(3) y(3) zd(3); 0 x(4) y(4) zd(4)]) / De;b(1)=det ([1 1 y(1) zd(1); 1 0 y(2) zd(2); 1 0 y(3) zd(3); 1 0 y(4) zd(4)]) / De; 
                  a(2)=det ([0 x(1) y(1) zd(1); 1 x(2) y(2) zd(2); 0 x(3) y(3) zd(3); 0 x(4) y(4) zd(4)]) / De;b(2)=det ([1 0 y(1) zd(1); 1 1 y(2) zd(2); 1 0 y(3) zd(3); 1 0 y(4) zd(4)]) / De; 
                  a(3)=det ([0 x(1) y(1) zd(1); 0 x(2) y(2) zd(2); 1 x(3) y(3) zd(3); 0 x(4) y(4) zd(4)]) / De;b(3)=det ([1 0 y(1) zd(1); 1 0 y(2) zd(2); 1 1 y(3) zd(3); 1 0 y(4) zd(4)]) / De;
                  a(4)=det ([0 x(1) y(1) zd(1); 0 x(2) y(2) zd(2); 0 x(3) y(3) zd(3); 1 x(4) y(4) zd(4)]) / De;b(4)=det ([1 0 y(1) zd(1); 1 0 y(2) zd(2); 1 0 y(3) zd(3); 1 1 y(4) zd(4)]) / De; 
                  c(1)=det ([1 x(1) 1 zd(1); 1 x(2) 0 zd(2); 1 x(3) 0 zd(3); 1 x(4) 0 zd(4)]) / De;d(1)=det ([1 x(1) y(1) 1; 1 x(2) y(2) 0; 1 x(3) y(3) 0; 1 x(4) y(4) 0]) / De;         
                  c(2)=det ([1 x(1) 0 zd(1); 1 x(2) 1 zd(2); 1 x(3) 0 zd(3); 1 x(4) 0 zd(4)]) / De;d(2)=det ([1 x(1) y(1) 0; 1 x(2) y(2) 1; 1 x(3) y(3) 0; 1 x(4) y(4) 0]) / De;         
                  c(3)=det ([1 x(1) 0 zd(1); 1 x(2) 0 zd(2); 1 x(3) 1 zd(3); 1 x(4) 0 zd(4)]) / De;d(3)=det ([1 x(1) y(1) 0; 1 x(2) y(2) 0; 1 x(3) y(3) 1; 1 x(4) y(4) 0]) / De;
                  c(4)=det ([1 x(1) 0 zd(1); 1 x(2) 0 zd(2); 1 x(3) 0 zd(3); 1 x(4) 1 zd(4)]) / De;d(4)=det ([1 x(1) y(1) 0; 1 x(2) y(2) 0; 1 x(3) y(3) 0; 1 x(4) y(4) 1]) / De;
                  z(1)=a(1)+b(1)*xx(ii)+c(1)*yy(jj)+d(1)*zz;z(2)=a(2)+b(2)*xx(ii)+c(2)*yy(jj)+d(2)*zz;z(3)=a(3)+b(3)*xx(ii)+c(3)*yy(jj)+d(3)*zz;z(4)=a(4)+b(4)*xx(ii)+c(4)*yy(jj)+d(4)*zz;
                  if (z(1)>0 || abs(z(1))<=error) && (z(2)>0 || abs(z(2))<=error) && (z(3)>0 || abs(z(3))<=error)&& (z(4)>0 || abs(z(4))<=error),ie=Is(aux);end
            end
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
subplot(2,3,1);pcolor(real(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,2);pcolor(real(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,3);pcolor(real(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,4);pcolor(real(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux');colorbar;axis tight;axis equal;
subplot(2,3,5);pcolor(real(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux');colorbar;axis tight;axis equal;
subplot(2,3,6);pcolor(real(Bz));shading interp;colormap(cmap);title('z Component Magnetic Flux');colorbar;axis tight;axis equal;



 
end




%{
figure
subplot(2,3,1);pcolor(real(Ex1));shading interp;hold on;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,2);pcolor(real(Ey1));shading interp;hold on;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,3);pcolor(real(Ez1));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,4);pcolor(real(Bx1));shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux');colorbar;axis tight;axis equal;
subplot(2,3,5);pcolor(real(By1));shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux');colorbar;axis tight;axis equal;
subplot(2,3,6);pcolor(real(Bz1));shading interp;colormap(cmap);title('z Component Magnetic Flux');colorbar;axis tight;axis equal;
%}
    %{
[cmap]=buildcmap('cbkry');figure;
subplot(2,3,1);surf(L1,L2,real(Ex'));shading interp;hold on;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,2);surf(L1,L2,real(Ey'));shading interp;hold on;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,3);surf(L1,L2,real(Ez'));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,4);surf(L1,L2,real(Bx'));shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux');colorbar;axis tight;%axis equal;
subplot(2,3,5);surf(L1,L2,real(By'));shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux');colorbar;axis tight;%axis equal;
subplot(2,3,6);surf(L1,L2,real(Bz'));shading interp;colormap(cmap);title('z Component Magnetic Flux');colorbar;axis tight;%axis equal;
    %}

%{
%{
figure;
Ex2=imagesc(real(Ex), 'Interpolation', 'bilinear');Ey2=imagesc(real(Ey), 'Interpolation', 'bilinear');Ez2=imagesc(real(Ez), 'Interpolation', 'bilinear');
Bx2=imagesc(real(Bx), 'Interpolation', 'bilinear');By2=imagesc(real(By), 'Interpolation', 'bilinear');Bz2=imagesc(real(Bz), 'Interpolation', 'bilinear');
subplot(2,3,1);pcolor(real(Ex2));shading interp;hold on;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,2);pcolor(real(Ey2));shading interp;hold on;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,3);pcolor(real(Ez2));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,4);pcolor(real(Bx2));shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux');colorbar;axis tight;axis equal;
subplot(2,3,5);pcolor(real(By2));shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux');colorbar;axis tight;axis equal;
subplot(2,3,6);pcolor(real(Bz2));shading interp;colormap(cmap);title('z Component Magnetic Flux');colorbar;axis tight;axis equal;
%}
%--------------------------------------- Smoothing ------------------------
%Ex=smoothdata(Ex);Ey=smoothdata(Ey);Ez=smoothdata(Ez);
%Ex=smoothdata(Ex,'rloess');Ey=smoothdata(Ey,'rloess');Ez=smoothdata(Ez,'rloess');
%Bx=smoothdata(Bx,'rloess');By=smoothdata(By,'rloess');Bz=smoothdata(Bz,'rloess');
Ex=smoothdata(Ex,'sgolay');Ey=smoothdata(Ey,'sgolay');Ez=smoothdata(Ez,'sgolay');
Bx=smoothdata(Bx,'sgolay');By=smoothdata(By,'sgolay');Bz=smoothdata(Bz,'sgolay');figure;
%--------------------------------------------------------------------------
subplot(2,3,1);pcolor(real(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,2);pcolor(real(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,3);pcolor(real(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,4);pcolor(real(Bx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux');colorbar;axis tight;axis equal;
subplot(2,3,5);pcolor(real(By));shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux');colorbar;axis tight;axis equal;
subplot(2,3,6);pcolor(real(Bz));shading interp;colormap(cmap);title('z Component Magnetic Flux');colorbar;axis tight;axis equal;
%--------------------------------------------------------------------------
figure;
En=sqrt(norm(Ex)^2 + norm(Ey)^2+norm(Ez)^2);Bn=sqrt(norm(Bx)^2 + norm(By)^2 + norm(Bz));Ex1=Ex/En;Ey1=Ey/En;Ez1=Ez/En;Bx1=Bx/Bn;By1=By/Bn;Bz1=Bz/En;
%K1=abs(Ex1)<1e-4;K2=abs(Ey1)<1e-4;K3=abs(Ez1)<1e-4;K4=abs(Bx1)<1e-6;K5=abs(By1)<1e-6;K6=abs(Bz1)<1e-6;
%Ex1(K1)=0;Ey1(K2)=0;Ez1(K3)=0;Bx1(K4)=0;By1(K5)=0;Bz1(K6)=0;
subplot(2,3,1);pcolor(real(Ex1));shading interp;hold on;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,2);pcolor(real(Ey1));shading interp;hold on;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,3);pcolor(real(Ez1));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
subplot(2,3,4);pcolor(real(Bx1));shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux');colorbar;axis tight;axis equal;
subplot(2,3,5);pcolor(real(By1));shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux');colorbar;axis tight;axis equal;
subplot(2,3,6);pcolor(real(Bz1));shading interp;colormap(cmap);title('z Component Magnetic Flux');colorbar;axis tight;axis equal;
%}
%}