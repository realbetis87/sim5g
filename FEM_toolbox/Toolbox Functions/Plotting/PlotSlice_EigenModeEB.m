%--------------------------------------------------------------------------
%{
                Plot E-B EigenMode Formulation EigenVectors

       1. X          : Uknown Field to be plotted (1 x Total Number Of  Uknowns)
       2. Position   : Plane Coordinates - for (-Ax) plane
                               Position=[-Ax;Nan;Nan;]
       3. Points     : [NumberOfPoints 1 Axis, NumberOfPOints 2 Axis]
       4. Axis       : Propagation Axis ("x","y","z")
       5. EigenValue : Complex Wavevector k
%}
%--------------------------------------------------------------------------
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz] = PlotSlice_EigenModeEB(tmodel,X,Position,Points,Axis,EigenValue),ElectromagneticConstants;
                  NumberOfFacets=numel(tmodel.Facets);NumberOfEdges=numel(tmodel.Edges);FacetVector=zeros(NumberOfFacets,1);EdgeVector=zeros(NumberOfEdges,1);
                  for ii=1:NumberOfFacets,facet=tmodel.Facets(ii);if(facet.UknownIndex~=0),FacetVector(ii)=-1i*X(abs(facet.UknownIndex))/c0;end,end
                  for ii=1:NumberOfEdges,edge=tmodel.Edges(ii);if(edge.UknownIndex~=0),EdgeVector(ii)=X(abs(edge.UknownIndex));end,end
                  Plane=find(~isnan(Position));N1=Points(1);N2=Points(2);empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
                  %Rectangular Domains
                  vertices=[tmodel.Vertices];xs=[vertices.X];ys=[vertices.Y];zs=[vertices.Z];
                  Xmin=min(xs);Ymin=min(ys);Zmin=min(zs);Xmax=max(xs);Ymax=max(ys);Zmax=max(zs);Wx=Xmax-Xmin;Wy=Ymax-Ymin;Wz=Zmax-Zmin;Wx=0;Wy=0;Wz=0;
                  Xmin=Xmin+Wx*1e-4;Xmax=Xmax-Wx/100;Ymin=Ymin+Wy*1e-4;Ymax=Ymax-Wy/100;Zmin=Zmin+Wz/100;Zmax=Zmax-1e-4;
                  switch Plane
                      case 1,xx=Position(1);yy=linspace(Ymin,Ymax,N1);zz=linspace(Zmin,Zmax,N2);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotXPlane(tmodel,EdgeVector,FacetVector,xx,yy,zz,N1,N2,Axis,EigenValue);
                      case 2,xx=linspace(Xmin,Xmax,N1);yy=Position(2);zz=linspace(Zmin,Zmax,N2);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotYPlane(tmodel,EdgeVector,FacetVector,xx,yy,zz,N1,N2,Axis,EigenValue);
                      case 3,xx=linspace(Xmin,Xmax,N1);yy=linspace(Ymin,Ymax,N2);zz=Position(3);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotZPlane(tmodel,EdgeVector,FacetVector,xx,yy,zz,N1,N2,Axis,EigenValue);
                   end
end
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotXPlane(tmodel,EdgeVector,FacetVector,xx,yy,zz,N1,N2,Axis,EigenValue),error=1000*eps;empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    Elements=[tmodel.Elements];Edges=[tmodel.Edges];Facets=[tmodel.Facets];Barycenter=[Elements.Barycenter];Ys=zeros(N1,N2);Zs=zeros(N1,N2);
    for ii=1:N1
        for jj=1:N2,distance=((xx-Barycenter(1,:)).^2+(yy(ii)-Barycenter(2,:)).^2 + (zz(jj)-Barycenter(3,:)).^2);[~,Is]=sort(distance);ie=0;aux=0;
            while(ie==0),aux=aux+1;element=Elements(Is(aux));a=element.As;b=element.Bs;c=element.Cs;d=element.Ds;
                z(1)=a(1)+b(1)*xx+c(1)*yy(ii)+d(1)*zz(jj);z(2)=a(2)+b(2)*xx+c(2)*yy(ii)+d(2)*zz(jj);z(3)=a(3)+b(3)*xx+c(3)*yy(ii)+d(3)*zz(jj);z(4)=a(4)+b(4)*xx+c(4)*yy(ii)+d(4)*zz(jj);
                if (z(1)>0 || abs(z(1))<=error) && (z(2)>0 || abs(z(2))<=error) && (z(3)>0 || abs(z(3))<=error)&& (z(4)>0 || abs(z(4))<=error),ie=Is(aux);end
            end,edges=[element.Edges];facets=[element.Facets];ll=[Edges(edges).Length];ss=[Facets(facets).Surface];Ys(ii,jj)=yy(ii);Zs(ii,jj)=zz(jj);
            wfx(1)=2*z(3)*(c(1)*d(2)-c(2)*d(1))+2*z(1)*(c(2)*d(3)-c(3)*d(2))+2*z(2)*(c(3)*d(1)-c(1)*d(3));wfy(1)=2*z(3)*(d(1)*b(2)-d(2)*b(1))+2*z(1)*(d(2)*b(3)-d(3)*b(2))+2*z(2)*(d(3)*b(1)-d(1)*b(3));wfz(1)=2*z(3)*(b(1)*c(2)-b(2)*c(1))+2*z(1)*(b(2)*c(3)-b(3)*c(2))+2*z(2)*(b(3)*c(1)-b(1)*c(3));
            wfx(2)=2*z(3)*(c(2)*d(4)-c(4)*d(2))+2*z(2)*(c(4)*d(3)-c(3)*d(4))+2*z(4)*(c(3)*d(2)-c(2)*d(3));wfy(2)=2*z(3)*(d(2)*b(4)-d(4)*b(2))+2*z(2)*(d(4)*b(3)-d(3)*b(4))+2*z(4)*(d(3)*b(2)-d(2)*b(3));wfz(2)=2*z(3)*(b(2)*c(4)-b(4)*c(2))+2*z(2)*(b(4)*c(3)-b(3)*c(4))+2*z(4)*(b(3)*c(2)-b(2)*c(3));
            wfx(3)=2*z(1)*(c(3)*d(4)-c(4)*d(3))+2*z(2)*(c(4)*d(1)-c(1)*d(4))+2*z(4)*(c(1)*d(3)-c(3)*d(1));wfy(3)=2*z(1)*(d(3)*b(4)-d(4)*b(3))+2*z(2)*(d(4)*b(1)-d(1)*b(4))+2*z(4)*(d(1)*b(3)-d(3)*b(1));wfz(3)=2*z(1)*(b(3)*c(4)-b(4)*c(3))+2*z(2)*(b(4)*c(1)-b(1)*c(4))+2*z(4)*(b(1)*c(3)-b(3)*c(1));
            wfx(4)=2*z(1)*(c(4)*d(2)-c(2)*d(4))+2*z(4)*(c(2)*d(1)-c(1)*d(2))+2*z(2)*(c(1)*d(4)-c(4)*d(1));wfy(4)=2*z(1)*(d(4)*b(2)-d(2)*b(4))+2*z(4)*(d(2)*b(1)-d(1)*b(2))+2*z(2)*(d(1)*b(4)-d(4)*b(1));wfz(4)=2*z(1)*(b(4)*c(2)-b(2)*c(4))+2*z(4)*(b(2)*c(1)-b(1)*c(2))+2*z(2)*(b(1)*c(4)-b(4)*c(1));
            wx(1)=z(1)*b(2)-z(2)*b(1); wx(2)=z(1)*b(3)-z(3)*b(1);wx(3)=z(1)*b(4)-z(4)*b(1);wx(4)=z(2)*b(3)-z(3)*b(2);wx(5)=z(2)*b(4)-z(4)*b(2);wx(6)=z(3)*b(4)-z(4)*b(3);wy(1)=z(1)*c(2)-z(2)*c(1); wy(2)=z(1)*c(3)-z(3)*c(1);wy(3)=z(1)*c(4)-z(4)*c(1);wy(4)=z(2)*c(3)-z(3)*c(2);wy(5)=z(2)*c(4)-z(4)*c(2);wy(6)=z(3)*c(4)-z(4)*c(3);
            wz(1)=z(1)*d(2)-z(2)*d(1); wz(2)=z(1)*d(3)-z(3)*d(1);wz(3)=z(1)*d(4)-z(4)*d(1);wz(4)=z(2)*d(3)-z(3)*d(2);wz(5)=z(2)*d(4)-z(4)*d(2);wz(6)=z(3)*d(4)-z(4)*d(3);
            %wx=wx.*ll;wy=wy.*ll;wz=wz.*ll;wfx=wfx.*ss;wfy=wfy.*ss;wfz=wfz.*ss;
            xp=xx;yp=yy(ii);zp=zz(jj);dfi=zeros(4,1);dei=zeros(6,1);oedges=[Edges(edges)];ofacets=[Facets(facets)];Ys(ii,jj)=yy(ii);Zs(ii,jj)=zz(jj);
            for mm=1:6,if(oedges(mm).UknownIndex==0),dei(mm)=1;elseif(sign(oedges(mm).UknownIndex)<0),dei(mm)=-1;else,dei(mm)=1;end,end
            for mm=1:4,if(ofacets(mm).UknownIndex==0),dfi(mm)=1;elseif(sign(ofacets(mm).UknownIndex)<0),dfi(mm)=-1;else,dfi(mm)=1;end,end
            switch Axis
                case "x"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*xp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*xp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*xp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);
                                                Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*xp);
                                                Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*xp);
                                                Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*xp);
                    end
                case "y"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*yp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*yp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*yp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*yp);
                                                         Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*yp);
                                                         Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*yp);
                    end
                case "z"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*sgn*dfi(kk)*FacetVector(facets(kk))*exp(-1i*EigenValue*zp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*sgn*dfi(kk)*FacetVector(facets(kk))*exp(-1i*EigenValue*zp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*sgn*dfi(kk)*FacetVector(facets(kk))*exp(-1i*EigenValue*zp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*zp);
                                                         Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*zp);
                                                         Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*zp);
                    end
            end
         end
    end
    rEx=real(Ex);[ii,jj]=find(abs(rEx)<1000*eps);rEx(ii,jj)=0;rEy=real(Ey);[ii,jj]=find(abs(rEy)<1000*eps);rEy(ii,jj)=0;rEz=real(Ez);[ii,jj]=find(abs(rEz)<1000*eps);rEz(ii,jj)=0;
    iEx=imag(Ex);[ii,jj]=find(abs(iEx)<1000*eps);iEx(ii,jj)=0;iEy=imag(Ey);[ii,jj]=find(abs(iEy)<1000*eps);iEy(ii,jj)=0;iEz=imag(Ez);[ii,jj]=find(abs(iEz)<1000*eps);iEz(ii,jj)=0;
    rBx=real(Bx);[ii,jj]=find(abs(rBx)<1000*eps);rBx(ii,jj)=0;rBy=real(By);[ii,jj]=find(abs(rBy)<1000*eps);rBy(ii,jj)=0;rBz=real(Bz);[ii,jj]=find(abs(rBz)<1000*eps);rBz(ii,jj)=0;
    iBx=imag(Bx);[ii,jj]=find(abs(iBx)<1000*eps);iBx(ii,jj)=0;iBy=imag(By);[ii,jj]=find(abs(iBy)<1000*eps);iBy(ii,jj)=0;iBz=imag(Bz);[ii,jj]=find(abs(iBz)<1000*eps);iBz(ii,jj)=0;
    [cmap]=buildcmap('cbkry');figure;
    subplot(1,3,1);pcolor(Ys,Zs,rEx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(Ys,Zs,rEy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(Ys,Zs,rEz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Electric Field Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(Ys,Zs,rBx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(Ys,Zs,rBy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(Ys,Zs,rBz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(Ys,Zs,iEx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(Ys,Zs,iEy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(Ys,Zs,iEz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(Ys,Zs,iBx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(Ys,Zs,iBy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(Ys,Zs,iBz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
end
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotYPlane(tmodel,EdgeVector,FacetVector,xx,yy,zz,N1,N2,Axis,EigenValue),error=1000*eps;empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    Elements=[tmodel.Elements];Edges=[tmodel.Edges];Facets=[tmodel.Facets];Barycenter=[Elements.Barycenter];Xs=zeros(N1,N2);Zs=zeros(N1,N2);
    for ii=1:N1
        for jj=1:N2,distance=((xx(ii)-Barycenter(1,:)).^2+(yy-Barycenter(2,:)).^2 + (zz(jj)-Barycenter(3,:)).^2);[~,Is]=sort(distance);ie=0;aux=0;
            while(ie==0),aux=aux+1;element=Elements(Is(aux));a=element.As;b=element.Bs;c=element.Cs;d=element.Ds;
                z(1)=a(1)+b(1)*xx(ii)+c(1)*yy+d(1)*zz(jj);z(2)=a(2)+b(2)*xx(ii)+c(2)*yy+d(2)*zz(jj);z(3)=a(3)+b(3)*xx(ii)+c(3)*yy+d(3)*zz(jj);z(4)=a(4)+b(4)*xx(ii)+c(4)*yy+d(4)*zz(jj);
                if (z(1)>0 || abs(z(1))<=error) && (z(2)>0 || abs(z(2))<=error) && (z(3)>0 || abs(z(3))<=error)&& (z(4)>0 || abs(z(4))<=error),ie=Is(aux);end
            end,edges=[element.Edges];facets=[element.Facets];ll=[Edges(edges).Length];ss=[Facets(facets).Surface];
            Xs(ii,jj)=xx(ii);Zs(ii,jj)=zz(jj);
            wfx(1)=2*z(3)*(c(1)*d(2)-c(2)*d(1))+2*z(1)*(c(2)*d(3)-c(3)*d(2))+2*z(2)*(c(3)*d(1)-c(1)*d(3));wfy(1)=2*z(3)*(d(1)*b(2)-d(2)*b(1))+2*z(1)*(d(2)*b(3)-d(3)*b(2))+2*z(2)*(d(3)*b(1)-d(1)*b(3));wfz(1)=2*z(3)*(b(1)*c(2)-b(2)*c(1))+2*z(1)*(b(2)*c(3)-b(3)*c(2))+2*z(2)*(b(3)*c(1)-b(1)*c(3));
            wfx(2)=2*z(3)*(c(2)*d(4)-c(4)*d(2))+2*z(2)*(c(4)*d(3)-c(3)*d(4))+2*z(4)*(c(3)*d(2)-c(2)*d(3));wfy(2)=2*z(3)*(d(2)*b(4)-d(4)*b(2))+2*z(2)*(d(4)*b(3)-d(3)*b(4))+2*z(4)*(d(3)*b(2)-d(2)*b(3));wfz(2)=2*z(3)*(b(2)*c(4)-b(4)*c(2))+2*z(2)*(b(4)*c(3)-b(3)*c(4))+2*z(4)*(b(3)*c(2)-b(2)*c(3));
            wfx(3)=2*z(1)*(c(3)*d(4)-c(4)*d(3))+2*z(2)*(c(4)*d(1)-c(1)*d(4))+2*z(4)*(c(1)*d(3)-c(3)*d(1));wfy(3)=2*z(1)*(d(3)*b(4)-d(4)*b(3))+2*z(2)*(d(4)*b(1)-d(1)*b(4))+2*z(4)*(d(1)*b(3)-d(3)*b(1));wfz(3)=2*z(1)*(b(3)*c(4)-b(4)*c(3))+2*z(2)*(b(4)*c(1)-b(1)*c(4))+2*z(4)*(b(1)*c(3)-b(3)*c(1));
            wfx(4)=2*z(1)*(c(4)*d(2)-c(2)*d(4))+2*z(4)*(c(2)*d(1)-c(1)*d(2))+2*z(2)*(c(1)*d(4)-c(4)*d(1));wfy(4)=2*z(1)*(d(4)*b(2)-d(2)*b(4))+2*z(4)*(d(2)*b(1)-d(1)*b(2))+2*z(2)*(d(1)*b(4)-d(4)*b(1));wfz(4)=2*z(1)*(b(4)*c(2)-b(2)*c(4))+2*z(4)*(b(2)*c(1)-b(1)*c(2))+2*z(2)*(b(1)*c(4)-b(4)*c(1));
            wx(1)=z(1)*b(2)-z(2)*b(1); wx(2)=z(1)*b(3)-z(3)*b(1);wx(3)=z(1)*b(4)-z(4)*b(1);wx(4)=z(2)*b(3)-z(3)*b(2);wx(5)=z(2)*b(4)-z(4)*b(2);wx(6)=z(3)*b(4)-z(4)*b(3);wy(1)=z(1)*c(2)-z(2)*c(1); wy(2)=z(1)*c(3)-z(3)*c(1);wy(3)=z(1)*c(4)-z(4)*c(1);wy(4)=z(2)*c(3)-z(3)*c(2);wy(5)=z(2)*c(4)-z(4)*c(2);wy(6)=z(3)*c(4)-z(4)*c(3);
            wz(1)=z(1)*d(2)-z(2)*d(1); wz(2)=z(1)*d(3)-z(3)*d(1);wz(3)=z(1)*d(4)-z(4)*d(1);wz(4)=z(2)*d(3)-z(3)*d(2);wz(5)=z(2)*d(4)-z(4)*d(2);wz(6)=z(3)*d(4)-z(4)*d(3);
            %wx=wx.*ll;wy=wy.*ll;wz=wz.*ll;wfx=wfx.*ss;wfy=wfy.*ss;wfz=wfz.*ss;
            xp=xx(ii);yp=yy;zp=zz(jj);dfi=zeros(4,1);dei=zeros(6,1);oedges=[Edges(edges)];ofacets=[Facets(facets)];
            for mm=1:6,if(oedges(mm).UknownIndex==0),dei(mm)=1;elseif(sign(oedges(mm).UknownIndex)<0),dei(mm)=-1;else,dei(mm)=1;end,end
            for mm=1:4,if(ofacets(mm).UknownIndex==0),dfi(mm)=1;elseif(sign(ofacets(mm).UknownIndex)<0),dfi(mm)=-1;else,dfi(mm)=1;end,end
            switch Axis
                case "x"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*xp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*xp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*xp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);
                                                Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*xp);
                                                Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*xp);
                                                Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*xp);
                    end
                case "y"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*yp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*yp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*yp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*yp);
                                                         Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*yp);
                                                         Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*yp);
                    end
                case "z"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*sgn*dfi(kk)*FacetVector(facets(kk))*exp(-1i*EigenValue*zp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*sgn*dfi(kk)*FacetVector(facets(kk))*exp(-1i*EigenValue*zp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*sgn*dfi(kk)*FacetVector(facets(kk))*exp(-1i*EigenValue*zp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*zp);
                                                         Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*zp);
                                                         Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*zp);
                    end
            end
        end
     end
    rEx=real(Ex);[ii,jj]=find(abs(rEx)<1000*eps);rEx(ii,jj)=0;rEy=real(Ey);[ii,jj]=find(abs(rEy)<1000*eps);rEy(ii,jj)=0;rEz=real(Ez);[ii,jj]=find(abs(rEz)<1000*eps);rEz(ii,jj)=0;
    iEx=imag(Ex);[ii,jj]=find(abs(iEx)<1000*eps);iEx(ii,jj)=0;iEy=imag(Ey);[ii,jj]=find(abs(iEy)<1000*eps);iEy(ii,jj)=0;iEz=imag(Ez);[ii,jj]=find(abs(iEz)<1000*eps);iEz(ii,jj)=0;
    rBx=real(Bx);[ii,jj]=find(abs(rBx)<1000*eps);rBx(ii,jj)=0;rBy=real(By);[ii,jj]=find(abs(rBy)<1000*eps);rBy(ii,jj)=0;rBz=real(Bz);[ii,jj]=find(abs(rBz)<1000*eps);rBz(ii,jj)=0;
    iBx=imag(Bx);[ii,jj]=find(abs(iBx)<1000*eps);iBx(ii,jj)=0;iBy=imag(By);[ii,jj]=find(abs(iBy)<1000*eps);iBy(ii,jj)=0;iBz=imag(Bz);[ii,jj]=find(abs(iBz)<1000*eps);iBz(ii,jj)=0;
    [cmap]=buildcmap('cbkry');figure;
    subplot(1,3,1);pcolor(Xs,Zs,rEx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(Xs,Zs,rEy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(Xs,Zs,rEz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Electric Field Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(Xs,Zs,rBx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(Xs,Zs,rBy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(Xs,Zs,rBz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(Xs,Zs,iEx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(Xs,Zs,iEy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(Xs,Zs,iEz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(Xs,Zs,iBx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(Xs,Zs,iBy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(Xs,Zs,iBz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
end
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotZPlane(tmodel,EdgeVector,FacetVector,xx,yy,zz,N1,N2,Axis,EigenValue),error=1000*eps;empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    Elements=[tmodel.Elements];Edges=[tmodel.Edges];Facets=[tmodel.Facets];Barycenter=[Elements.Barycenter];Xs=zeros(N1,N2);Ys=zeros(N1,N2);
    for ii=1:N1
        for jj=1:N2,distance=((xx(ii)-Barycenter(1,:)).^2+(yy(jj)-Barycenter(2,:)).^2 + (zz-Barycenter(3,:)).^2);[~,Is]=sort(distance);ie=0;aux=0;
            while(ie==0),aux=aux+1;element=Elements(Is(aux));a=element.As;b=element.Bs;c=element.Cs;d=element.Ds;
                z(1)=a(1)+b(1)*xx(ii)+c(1)*yy(jj)+d(1)*zz;
                z(2)=a(2)+b(2)*xx(ii)+c(2)*yy(jj)+d(2)*zz;
                z(3)=a(3)+b(3)*xx(ii)+c(3)*yy(jj)+d(3)*zz;
                z(4)=a(4)+b(4)*xx(ii)+c(4)*yy(jj)+d(4)*zz;
                if (z(1)>0 || abs(z(1))<=error) && (z(2)>0 || abs(z(2))<=error) && (z(3)>0 || abs(z(3))<=error)&& (z(4)>0 || abs(z(4))<=error),ie=Is(aux);end
            end,edges=[element.Edges];facets=[element.Facets];
            Xs(ii,jj)=xx(ii);Ys(ii,jj)=yy(jj);
            wfx(1)=2*z(3)*(c(1)*d(2)-c(2)*d(1))+2*z(1)*(c(2)*d(3)-c(3)*d(2))+2*z(2)*(c(3)*d(1)-c(1)*d(3));wfy(1)=2*z(3)*(d(1)*b(2)-d(2)*b(1))+2*z(1)*(d(2)*b(3)-d(3)*b(2))+2*z(2)*(d(3)*b(1)-d(1)*b(3));wfz(1)=2*z(3)*(b(1)*c(2)-b(2)*c(1))+2*z(1)*(b(2)*c(3)-b(3)*c(2))+2*z(2)*(b(3)*c(1)-b(1)*c(3));
            wfx(2)=2*z(3)*(c(2)*d(4)-c(4)*d(2))+2*z(2)*(c(4)*d(3)-c(3)*d(4))+2*z(4)*(c(3)*d(2)-c(2)*d(3));wfy(2)=2*z(3)*(d(2)*b(4)-d(4)*b(2))+2*z(2)*(d(4)*b(3)-d(3)*b(4))+2*z(4)*(d(3)*b(2)-d(2)*b(3));wfz(2)=2*z(3)*(b(2)*c(4)-b(4)*c(2))+2*z(2)*(b(4)*c(3)-b(3)*c(4))+2*z(4)*(b(3)*c(2)-b(2)*c(3));
            wfx(3)=2*z(1)*(c(3)*d(4)-c(4)*d(3))+2*z(2)*(c(4)*d(1)-c(1)*d(4))+2*z(4)*(c(1)*d(3)-c(3)*d(1));wfy(3)=2*z(1)*(d(3)*b(4)-d(4)*b(3))+2*z(2)*(d(4)*b(1)-d(1)*b(4))+2*z(4)*(d(1)*b(3)-d(3)*b(1));wfz(3)=2*z(1)*(b(3)*c(4)-b(4)*c(3))+2*z(2)*(b(4)*c(1)-b(1)*c(4))+2*z(4)*(b(1)*c(3)-b(3)*c(1));
            wfx(4)=2*z(1)*(c(4)*d(2)-c(2)*d(4))+2*z(4)*(c(2)*d(1)-c(1)*d(2))+2*z(2)*(c(1)*d(4)-c(4)*d(1));wfy(4)=2*z(1)*(d(4)*b(2)-d(2)*b(4))+2*z(4)*(d(2)*b(1)-d(1)*b(2))+2*z(2)*(d(1)*b(4)-d(4)*b(1));wfz(4)=2*z(1)*(b(4)*c(2)-b(2)*c(4))+2*z(4)*(b(2)*c(1)-b(1)*c(2))+2*z(2)*(b(1)*c(4)-b(4)*c(1));
            wx(1)=z(1)*b(2)-z(2)*b(1); wx(2)=z(1)*b(3)-z(3)*b(1);wx(3)=z(1)*b(4)-z(4)*b(1);wx(4)=z(2)*b(3)-z(3)*b(2);wx(5)=z(2)*b(4)-z(4)*b(2);wx(6)=z(3)*b(4)-z(4)*b(3);wy(1)=z(1)*c(2)-z(2)*c(1); wy(2)=z(1)*c(3)-z(3)*c(1);wy(3)=z(1)*c(4)-z(4)*c(1);wy(4)=z(2)*c(3)-z(3)*c(2);wy(5)=z(2)*c(4)-z(4)*c(2);wy(6)=z(3)*c(4)-z(4)*c(3);
            wz(1)=z(1)*d(2)-z(2)*d(1); wz(2)=z(1)*d(3)-z(3)*d(1);wz(3)=z(1)*d(4)-z(4)*d(1);wz(4)=z(2)*d(3)-z(3)*d(2);wz(5)=z(2)*d(4)-z(4)*d(2);wz(6)=z(3)*d(4)-z(4)*d(3);
            
            xp=xx(ii);yp=yy(jj);zp=zz;dfi=zeros(4,1);dei=zeros(6,1);oedges=[Edges(edges)];ofacets=[Facets(facets)];edges=[element.Edges];facets=[element.Facets];
            ll=[oedges.Length];ss=[ofacets.Surface];ll=ones(6,1);ss=ones(4,1);
            wx=wx.*ll;wy=wy.*ll;wz=wz.*ll;wfx=wfx.*ss;wfy=wfy.*ss;wfz=wfz.*ss;%dsi=zeros(6,1);
            for mm=1:6,if(oedges(mm).UknownIndex==0),dei(mm)=1;elseif(sign(oedges(mm).UknownIndex)<0),dei(mm)=-1;else,dei(mm)=1;end,end
            for mm=1:4,if(ofacets(mm).UknownIndex==0),dfi(mm)=1;elseif(sign(ofacets(mm).UknownIndex)<0),dfi(mm)=-1;else,dfi(mm)=1;end,end
            switch Axis
                case "x"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*xp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*xp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*xp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);
                                                Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*xp);
                                                Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*xp);
                                                Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*xp);
                    end
                case "y"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*yp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*yp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk))*exp(-1i*EigenValue*yp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*yp);
                                                         Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*yp);
                                                         Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*yp);
                    end
                case "z"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*sgn*dfi(kk)*FacetVector(facets(kk))*exp(-1i*EigenValue*zp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*sgn*dfi(kk)*FacetVector(facets(kk))*exp(-1i*EigenValue*zp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*sgn*dfi(kk)*FacetVector(facets(kk))*exp(-1i*EigenValue*zp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*zp);
                                                         Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*zp);
                                                         Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk))*exp(-1i*EigenValue*zp);
                    end
            end
        end
    end
    rEx=real(Ex);[ii,jj]=find(abs(rEx)<1000*eps);rEx(ii,jj)=0;rEy=real(Ey);[ii,jj]=find(abs(rEy)<1000*eps);rEy(ii,jj)=0;rEz=real(Ez);[ii,jj]=find(abs(rEz)<1000*eps);rEz(ii,jj)=0;
    iEx=imag(Ex);[ii,jj]=find(abs(iEx)<1000*eps);iEx(ii,jj)=0;iEy=imag(Ey);[ii,jj]=find(abs(iEy)<1000*eps);iEy(ii,jj)=0;iEz=imag(Ez);[ii,jj]=find(abs(iEz)<1000*eps);iEz(ii,jj)=0;
    rBx=real(Bx);[ii,jj]=find(abs(rBx)<1000*eps);rBx(ii,jj)=0;rBy=real(By);[ii,jj]=find(abs(rBy)<1000*eps);rBy(ii,jj)=0;rBz=real(Bz);[ii,jj]=find(abs(rBz)<1000*eps);rBz(ii,jj)=0;
    iBx=imag(Bx);[ii,jj]=find(abs(iBx)<1000*eps);iBx(ii,jj)=0;iBy=imag(By);[ii,jj]=find(abs(iBy)<1000*eps);iBy(ii,jj)=0;iBz=imag(Bz);[ii,jj]=find(abs(iBz)<1000*eps);iBz(ii,jj)=0;
    [cmap]=buildcmap('cbkry');figure;
    subplot(1,3,1);pcolor(Xs,Ys,rEx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(Xs,Ys,rEy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(Xs,Ys,rEz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Electric Field Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(Xs,Ys,rBx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(Xs,Ys,rBy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(Xs,Ys,rBz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(Xs,Ys,iEx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(Xs,Ys,iEy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(Xs,Ys,iEz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(Xs,Ys,iBx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(Xs,Ys,iBy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(Xs,Ys,iBz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
end