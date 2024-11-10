%--------------------------------------------------------------------------
%{
            Plot E-B EigenMode Formulation Eigen Vectors 

       1. X          : Uknown Field to be plotted (1 x Total Number Of  Uknowns)
       2. Position   : Plane Coordinates - for (-Ax) plane
                               Position=[-Ax;Nan;Nan;]
       3. Points     : [NumberOfPoints 1 Axis, NumberOfPOints 2 Axis]
       4. Axis       : Propagation Axis ("x","y","z")
       5. EigenValue : Complex Wavevector k
%}
%--------------------------------------------------------------------------
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz] = PlotPlane_EB_EigenMode(TModel,X,Position,Points,Axis,EigenValue),ElectromagneticConstants;
    NumberOfFacets=numel(TModel.Facets);NumberOfEdges=numel(TModel.Edges);FacetVector=zeros(NumberOfFacets,1);EdgeVector=zeros(NumberOfEdges,1);
    for ii=1:NumberOfFacets,facet=TModel.Facets(ii);if(facet.UknownIndex~=0),FacetVector(ii)=X(abs(abs(facet.UknownIndex)));end,end%FacetVector(ii)=-1i*X(abs(facet.UknownIndex))/c0;end,end
    for ii=1:NumberOfEdges,edge=TModel.Edges(ii);if(edge.UknownIndex~=0),EdgeVector(ii)=X(abs(edge.UknownIndex));end,end
    Plane=find(~isnan(Position));N1=Points(1);N2=Points(2);empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    vertices=[TModel.Vertices];xs=[vertices.X];ys=[vertices.Y];zs=[vertices.Z];
    Xmin=min(xs);Ymin=min(ys);Zmin=min(zs);Xmax=max(xs);Ymax=max(ys);Zmax=max(zs);
    switch Plane
        case 1,xx=Position(1);y=linspace(Ymin-Ymin*1e-3,Ymax-Ymax*1e-3,N1);z=linspace(Zmin-Zmin*1e-3,Zmax-Zmax*1e-3,N2);[yy,zz]=meshgrid(y,z);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotXPlane(TModel,EdgeVector,FacetVector,xx,yy,zz,Axis,EigenValue);
        case 2,x=linspace(Xmin-Xmin*1e-3,Xmax-Xmax*1e-3,N1);yy=Position(2);z=linspace(Zmin-Zmin*1e-3,Zmax-Zmax*1e-3,N2);[xx,zz]=meshgrid(x,z);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotYPlane(TModel,EdgeVector,FacetVector,xx,yy,zz,Axis,EigenValue);
        case 3,x=linspace(Xmin-Xmin*1e-3,Xmax-Xmax*1e-3,N1);y=linspace(Ymin-Ymin*1e-3,Ymax-Ymax*1e-3,N2);zz=Position(3);[xx,yy]=meshgrid(x,y);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotZPlane(TModel,EdgeVector,FacetVector,xx,yy,zz,Axis,EigenValue);
     end
end
%--------------------------------------------------------------------------
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotXPlane(TModel,EdgeVector,FacetVector,xx,yy,zz,Axis,EigenValue),error=1000*eps;empty=zeros(size(yy,1),size(yy,2));Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    Elements=[TModel.Elements];Barycenter=[Elements.Barycenter];
    for ii=1:size(yy,1)
        for jj=1:size(yy,2),distance=((xx-Barycenter(1,:)).^2+(yy(ii,jj)-Barycenter(2,:)).^2 + (zz(ii,jj)-Barycenter(3,:)).^2);[~,Is]=sort(distance);ie=0;aux=0;
            while(ie==0),aux=aux+1;element=Elements(Is(aux));a=element.As;b=element.Bs;c=element.Cs;d=element.Ds;
                zeta(1)=a(1)+b(1)*xx+c(1)*yy(ii,jj)+d(1)*zz(ii,jj);
                zeta(2)=a(2)+b(2)*xx+c(2)*yy(ii,jj)+d(2)*zz(ii,jj);
                zeta(3)=a(3)+b(3)*xx+c(3)*yy(ii,jj)+d(3)*zz(ii,jj);
                zeta(4)=a(4)+b(4)*xx+c(4)*yy(ii,jj)+d(4)*zz(ii,jj);
                if (zeta(1)>0 || abs(zeta(1))<=error) && (zeta(2)>0 || abs(zeta(2))<=error) && (zeta(3)>0 || abs(zeta(3))<=error)&& (zeta(4)>0 || abs(zeta(4))<=error),ie=Is(aux);end
            end,edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];ll=[edges.Length];ss=[facets.Surface];

            wfx(1)=2*zeta(3)*(c(1)*d(2)-c(2)*d(1))+2*zeta(1)*(c(2)*d(3)-c(3)*d(2))+2*zeta(2)*(c(3)*d(1)-c(1)*d(3));
            wfx(2)=2*zeta(3)*(c(2)*d(4)-c(4)*d(2))+2*zeta(2)*(c(4)*d(3)-c(3)*d(4))+2*zeta(4)*(c(3)*d(2)-c(2)*d(3));
            wfx(3)=2*zeta(1)*(c(3)*d(4)-c(4)*d(3))+2*zeta(2)*(c(4)*d(1)-c(1)*d(4))+2*zeta(4)*(c(1)*d(3)-c(3)*d(1));
            wfx(4)=2*zeta(1)*(c(4)*d(2)-c(2)*d(4))+2*zeta(4)*(c(2)*d(1)-c(1)*d(2))+2*zeta(2)*(c(1)*d(4)-c(4)*d(1));
            
            wfy(1)=2*zeta(3)*(d(1)*b(2)-d(2)*b(1))+2*zeta(1)*(d(2)*b(3)-d(3)*b(2))+2*zeta(2)*(d(3)*b(1)-d(1)*b(3));
            wfy(2)=2*zeta(3)*(d(2)*b(4)-d(4)*b(2))+2*zeta(2)*(d(4)*b(3)-d(3)*b(4))+2*zeta(4)*(d(3)*b(2)-d(2)*b(3));
            wfy(3)=2*zeta(1)*(d(3)*b(4)-d(4)*b(3))+2*zeta(2)*(d(4)*b(1)-d(1)*b(4))+2*zeta(4)*(d(1)*b(3)-d(3)*b(1));
            wfy(4)=2*zeta(1)*(d(4)*b(2)-d(2)*b(4))+2*zeta(4)*(d(2)*b(1)-d(1)*b(2))+2*zeta(2)*(d(1)*b(4)-d(4)*b(1));

            wfz(1)=2*zeta(3)*(b(1)*c(2)-b(2)*c(1))+2*zeta(1)*(b(2)*c(3)-b(3)*c(2))+2*zeta(2)*(b(3)*c(1)-b(1)*c(3));
            wfz(2)=2*zeta(3)*(b(2)*c(4)-b(4)*c(2))+2*zeta(2)*(b(4)*c(3)-b(3)*c(4))+2*zeta(4)*(b(3)*c(2)-b(2)*c(3));
            wfz(3)=2*zeta(1)*(b(3)*c(4)-b(4)*c(3))+2*zeta(2)*(b(4)*c(1)-b(1)*c(4))+2*zeta(4)*(b(1)*c(3)-b(3)*c(1));
            wfz(4)=2*zeta(1)*(b(4)*c(2)-b(2)*c(4))+2*zeta(4)*(b(2)*c(1)-b(1)*c(2))+2*zeta(2)*(b(1)*c(4)-b(4)*c(1));

            wx(1)=zeta(1)*b(2)-zeta(2)*b(1);      wy(1)=zeta(1)*c(2)-zeta(2)*c(1);      wz(1)=zeta(1)*d(2)-zeta(2)*d(1); 
            wx(2)=zeta(1)*b(3)-zeta(3)*b(1);      wy(2)=zeta(1)*c(3)-zeta(3)*c(1);      wz(2)=zeta(1)*d(3)-zeta(3)*d(1);
            wx(3)=zeta(1)*b(4)-zeta(4)*b(1);      wy(3)=zeta(1)*c(4)-zeta(4)*c(1);      wz(3)=zeta(1)*d(4)-zeta(4)*d(1);
            wx(4)=zeta(2)*b(3)-zeta(3)*b(2);      wy(4)=zeta(2)*c(3)-zeta(3)*c(2);      wz(4)=zeta(2)*d(3)-zeta(3)*d(2);
            wx(5)=zeta(2)*b(4)-zeta(4)*b(2);      wy(5)=zeta(2)*c(4)-zeta(4)*c(2);      wz(5)=zeta(2)*d(4)-zeta(4)*d(2);
            wx(6)=zeta(3)*b(4)-zeta(4)*b(3);      wy(6)=zeta(3)*c(4)-zeta(4)*c(3);      wz(6)=zeta(3)*d(4)-zeta(4)*d(3);
             
            wx=wx.*ll;wy=wy.*ll;wz=wz.*ll;wfx=wfx.*ss;wfy=wfy.*ss;wfz=wfz.*ss;

            xp=xx;yp=yy(ii,jj);zp=zz(ii,jj);dfi=zeros(4,1);dei=zeros(6,1);
            for mm=1:6,if(edges(mm).UknownIndex==0),dei(mm)=1;elseif(sign(edges(mm).UknownIndex)<0),dei(mm)=-1;else,dei(mm)=1;end,end
            for mm=1:4,if(facets(mm).UknownIndex==0),dfi(mm)=1;elseif(sign(facets(mm).UknownIndex)<0),dfi(mm)=-1;else,dfi(mm)=1;end,end
            switch Axis
                case "x"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*xp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*xp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*xp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);
                                                Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*xp);
                                                Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*xp);
                                                Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*xp);
                    end
                case "y"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*yp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*yp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*yp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*yp);
                                                         Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*yp);
                                                         Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*yp);
                    end
                case "z"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*sgn*dfi(kk)*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*zp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*sgn*dfi(kk)*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*zp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*sgn*dfi(kk)*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*zp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*zp);
                                                         Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*zp);
                                                         Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*zp);
                    end
            end
         end
    end
    rEx=real(Ex);[ii,jj]=find(abs(rEx)<100*eps);rEx(ii,jj)=0;rEy=real(Ey);[ii,jj]=find(abs(rEy)<10*eps);rEy(ii,jj)=0;rEz=real(Ez);[ii,jj]=find(abs(rEz)<10*eps);rEz(ii,jj)=0;
    iEx=imag(Ex);[ii,jj]=find(abs(iEx)<10*eps);iEx(ii,jj)=0;iEy=imag(Ey);[ii,jj]=find(abs(iEy)<10*eps);iEy(ii,jj)=0;iEz=imag(Ez);[ii,jj]=find(abs(iEz)<10*eps);iEz(ii,jj)=0;
    rBx=real(Bx);[ii,jj]=find(abs(rBx)<10*eps);rBx(ii,jj)=0;rBy=real(By);[ii,jj]=find(abs(rBy)<10*eps);rBy(ii,jj)=0;rBz=real(Bz);[ii,jj]=find(abs(rBz)<10*eps);rBz(ii,jj)=0;
    iBx=imag(Bx);[ii,jj]=find(abs(iBx)<10*eps);iBx(ii,jj)=0;iBy=imag(By);[ii,jj]=find(abs(iBy)<10*eps);iBy(ii,jj)=0;iBz=imag(Bz);[ii,jj]=find(abs(iBz)<10*eps);iBz(ii,jj)=0;
    [cmap]=buildcmap('cbkry');figure;
    subplot(1,3,1);pcolor(yy,zz,rEx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(yy,zz,rEy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(yy,zz,rEz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Electric Field Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(yy,zz,rBx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(yy,zz,rBy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(yy,zz,rBz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(yy,zz,iEx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(yy,zz,iEy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(yy,zz,iEz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(yy,zz,iBx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(yy,zz,iBy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(yy,zz,iBz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
end
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotYPlane(TModel,EdgeVector,FacetVector,xx,yy,zz,Axis,EigenValue),error=1000*eps;empty=zeros(size(xx,1),size(xx,2));Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    Elements=[TModel.Elements];Barycenter=[Elements.Barycenter];
    for ii=1:size(xx,1)
        for jj=1:size(xx,2),distance=((xx(ii,jj)-Barycenter(1,:)).^2+(yy-Barycenter(2,:)).^2 + (zz(ii,jj)-Barycenter(3,:)).^2);[~,Is]=sort(distance);ie=0;aux=0;
            while(ie==0),aux=aux+1;element=Elements(Is(aux));a=element.As;b=element.Bs;c=element.Cs;d=element.Ds;
                zeta(1)=a(1)+b(1)*xx(ii,jj)+c(1)*yy+d(1)*zz(ii,jj);
                zeta(2)=a(2)+b(2)*xx(ii,jj)+c(2)*yy+d(2)*zz(ii,jj);
                zeta(3)=a(3)+b(3)*xx(ii,jj)+c(3)*yy+d(3)*zz(ii,jj);
                zeta(4)=a(4)+b(4)*xx(ii,jj)+c(4)*yy+d(4)*zz(ii,jj);
                if (zeta(1)>0 || abs(zeta(1))<=error) && (zeta(2)>0 || abs(zeta(2))<=error) && (zeta(3)>0 || abs(zeta(3))<=error)&& (zeta(4)>0 || abs(zeta(4))<=error),ie=Is(aux);end
            end,edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];ll=[edges.Length];ss=[facets.Surface];

            wfx(1)=2*zeta(3)*(c(1)*d(2)-c(2)*d(1))+2*zeta(1)*(c(2)*d(3)-c(3)*d(2))+2*zeta(2)*(c(3)*d(1)-c(1)*d(3));
            wfx(2)=2*zeta(3)*(c(2)*d(4)-c(4)*d(2))+2*zeta(2)*(c(4)*d(3)-c(3)*d(4))+2*zeta(4)*(c(3)*d(2)-c(2)*d(3));
            wfx(3)=2*zeta(1)*(c(3)*d(4)-c(4)*d(3))+2*zeta(2)*(c(4)*d(1)-c(1)*d(4))+2*zeta(4)*(c(1)*d(3)-c(3)*d(1));
            wfx(4)=2*zeta(1)*(c(4)*d(2)-c(2)*d(4))+2*zeta(4)*(c(2)*d(1)-c(1)*d(2))+2*zeta(2)*(c(1)*d(4)-c(4)*d(1));
            
            wfy(1)=2*zeta(3)*(d(1)*b(2)-d(2)*b(1))+2*zeta(1)*(d(2)*b(3)-d(3)*b(2))+2*zeta(2)*(d(3)*b(1)-d(1)*b(3));
            wfy(2)=2*zeta(3)*(d(2)*b(4)-d(4)*b(2))+2*zeta(2)*(d(4)*b(3)-d(3)*b(4))+2*zeta(4)*(d(3)*b(2)-d(2)*b(3));
            wfy(3)=2*zeta(1)*(d(3)*b(4)-d(4)*b(3))+2*zeta(2)*(d(4)*b(1)-d(1)*b(4))+2*zeta(4)*(d(1)*b(3)-d(3)*b(1));
            wfy(4)=2*zeta(1)*(d(4)*b(2)-d(2)*b(4))+2*zeta(4)*(d(2)*b(1)-d(1)*b(2))+2*zeta(2)*(d(1)*b(4)-d(4)*b(1));

            wfz(1)=2*zeta(3)*(b(1)*c(2)-b(2)*c(1))+2*zeta(1)*(b(2)*c(3)-b(3)*c(2))+2*zeta(2)*(b(3)*c(1)-b(1)*c(3));
            wfz(2)=2*zeta(3)*(b(2)*c(4)-b(4)*c(2))+2*zeta(2)*(b(4)*c(3)-b(3)*c(4))+2*zeta(4)*(b(3)*c(2)-b(2)*c(3));
            wfz(3)=2*zeta(1)*(b(3)*c(4)-b(4)*c(3))+2*zeta(2)*(b(4)*c(1)-b(1)*c(4))+2*zeta(4)*(b(1)*c(3)-b(3)*c(1));
            wfz(4)=2*zeta(1)*(b(4)*c(2)-b(2)*c(4))+2*zeta(4)*(b(2)*c(1)-b(1)*c(2))+2*zeta(2)*(b(1)*c(4)-b(4)*c(1));

            wx(1)=zeta(1)*b(2)-zeta(2)*b(1);      wy(1)=zeta(1)*c(2)-zeta(2)*c(1);      wz(1)=zeta(1)*d(2)-zeta(2)*d(1); 
            wx(2)=zeta(1)*b(3)-zeta(3)*b(1);      wy(2)=zeta(1)*c(3)-zeta(3)*c(1);      wz(2)=zeta(1)*d(3)-zeta(3)*d(1);
            wx(3)=zeta(1)*b(4)-zeta(4)*b(1);      wy(3)=zeta(1)*c(4)-zeta(4)*c(1);      wz(3)=zeta(1)*d(4)-zeta(4)*d(1);
            wx(4)=zeta(2)*b(3)-zeta(3)*b(2);      wy(4)=zeta(2)*c(3)-zeta(3)*c(2);      wz(4)=zeta(2)*d(3)-zeta(3)*d(2);
            wx(5)=zeta(2)*b(4)-zeta(4)*b(2);      wy(5)=zeta(2)*c(4)-zeta(4)*c(2);      wz(5)=zeta(2)*d(4)-zeta(4)*d(2);
            wx(6)=zeta(3)*b(4)-zeta(4)*b(3);      wy(6)=zeta(3)*c(4)-zeta(4)*c(3);      wz(6)=zeta(3)*d(4)-zeta(4)*d(3);
             
            wx=wx.*ll;wy=wy.*ll;wz=wz.*ll;wfx=wfx.*ss;wfy=wfy.*ss;wfz=wfz.*ss;

            xp=xx(ii,jj);yp=yy;zp=zz(ii,jj);dfi=zeros(4,1);dei=zeros(6,1);
            for mm=1:6,if(edges(mm).UknownIndex==0),dei(mm)=1;elseif(sign(edges(mm).UknownIndex)<0),dei(mm)=-1;else,dei(mm)=1;end,end
            for mm=1:4,if(facets(mm).UknownIndex==0),dfi(mm)=1;elseif(sign(facets(mm).UknownIndex)<0),dfi(mm)=-1;else,dfi(mm)=1;end,end
            switch Axis
                case "x"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*xp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*xp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*xp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);
                                                Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*xp);
                                                Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*xp);
                                                Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*xp);
                    end
                case "y"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*yp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*yp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*yp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*yp);
                                                         Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*yp);
                                                         Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*yp);
                    end
                case "z"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*sgn*dfi(kk)*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*zp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*sgn*dfi(kk)*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*zp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*sgn*dfi(kk)*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*zp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*zp);
                                                         Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*zp);
                                                         Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*zp);
                    end
            end
         end
    end
    rEx=real(Ex);[ii,jj]=find(abs(rEx)<1000*eps);rEx(ii,jj)=0;rEy=real(Ey);[ii,jj]=find(abs(rEy)<1000*eps);rEy(ii,jj)=0;rEz=real(Ez);[ii,jj]=find(abs(rEz)<1000*eps);rEz(ii,jj)=0;
    iEx=imag(Ex);[ii,jj]=find(abs(iEx)<1000*eps);iEx(ii,jj)=0;iEy=imag(Ey);[ii,jj]=find(abs(iEy)<1000*eps);iEy(ii,jj)=0;iEz=imag(Ez);[ii,jj]=find(abs(iEz)<1000*eps);iEz(ii,jj)=0;
    rBx=real(Bx);[ii,jj]=find(abs(rBx)<1000*eps);rBx(ii,jj)=0;rBy=real(By);[ii,jj]=find(abs(rBy)<1000*eps);rBy(ii,jj)=0;rBz=real(Bz);[ii,jj]=find(abs(rBz)<1000*eps);rBz(ii,jj)=0;
    iBx=imag(Bx);[ii,jj]=find(abs(iBx)<1000*eps);iBx(ii,jj)=0;iBy=imag(By);[ii,jj]=find(abs(iBy)<1000*eps);iBy(ii,jj)=0;iBz=imag(Bz);[ii,jj]=find(abs(iBz)<1000*eps);iBz(ii,jj)=0;
    [cmap]=buildcmap('cbkry');figure;
    subplot(1,3,1);pcolor(xx,zz,rEx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(xx,zz,rEy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(xx,zz,rEz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Electric Field Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(xx,zz,rBx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(xx,zz,rBy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(xx,zz,rBz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(xx,zz,iEx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(xx,zz,iEy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(xx,zz,iEz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(xx,zz,iBx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(xx,zz,iBy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(xx,zz,iBz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
end
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotZPlane(TModel,EdgeVector,FacetVector,xx,yy,zz,Axis,EigenValue),error=1000*eps;empty=zeros(size(xx,1),size(xx,2));Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    Elements=[TModel.Elements];Barycenter=[Elements.Barycenter];
    for ii=1:size(xx,1)
        for jj=1:size(xx,2),distance=((xx(ii,jj)-Barycenter(1,:)).^2+(yy(ii,jj)-Barycenter(2,:)).^2 + (zz-Barycenter(3,:)).^2);[~,Is]=sort(distance);ie=0;aux=0;
            while(ie==0),aux=aux+1;element=Elements(Is(aux));a=element.As;b=element.Bs;c=element.Cs;d=element.Ds;
                zeta(1)=a(1)+b(1)*xx(ii,jj)+c(1)*yy(ii,jj)+d(1)*zz;
                zeta(2)=a(2)+b(2)*xx(ii,jj)+c(2)*yy(ii,jj)+d(2)*zz;
                zeta(3)=a(3)+b(3)*xx(ii,jj)+c(3)*yy(ii,jj)+d(3)*zz;
                zeta(4)=a(4)+b(4)*xx(ii,jj)+c(4)*yy(ii,jj)+d(4)*zz;
                if (zeta(1)>0 || abs(zeta(1))<=error) && (zeta(2)>0 || abs(zeta(2))<=error) && (zeta(3)>0 || abs(zeta(3))<=error)&& (zeta(4)>0 || abs(zeta(4))<=error),ie=Is(aux);end
            end,edges=[TModel.Edges(element.Edges)];facets=[TModel.Facets(element.Facets)];ll=[edges.Length];ss=[facets.Surface];

            wfx(1)=2*zeta(3)*(c(1)*d(2)-c(2)*d(1))+2*zeta(1)*(c(2)*d(3)-c(3)*d(2))+2*zeta(2)*(c(3)*d(1)-c(1)*d(3));
            wfx(2)=2*zeta(3)*(c(2)*d(4)-c(4)*d(2))+2*zeta(2)*(c(4)*d(3)-c(3)*d(4))+2*zeta(4)*(c(3)*d(2)-c(2)*d(3));
            wfx(3)=2*zeta(1)*(c(3)*d(4)-c(4)*d(3))+2*zeta(2)*(c(4)*d(1)-c(1)*d(4))+2*zeta(4)*(c(1)*d(3)-c(3)*d(1));
            wfx(4)=2*zeta(1)*(c(4)*d(2)-c(2)*d(4))+2*zeta(4)*(c(2)*d(1)-c(1)*d(2))+2*zeta(2)*(c(1)*d(4)-c(4)*d(1));
            
            wfy(1)=2*zeta(3)*(d(1)*b(2)-d(2)*b(1))+2*zeta(1)*(d(2)*b(3)-d(3)*b(2))+2*zeta(2)*(d(3)*b(1)-d(1)*b(3));
            wfy(2)=2*zeta(3)*(d(2)*b(4)-d(4)*b(2))+2*zeta(2)*(d(4)*b(3)-d(3)*b(4))+2*zeta(4)*(d(3)*b(2)-d(2)*b(3));
            wfy(3)=2*zeta(1)*(d(3)*b(4)-d(4)*b(3))+2*zeta(2)*(d(4)*b(1)-d(1)*b(4))+2*zeta(4)*(d(1)*b(3)-d(3)*b(1));
            wfy(4)=2*zeta(1)*(d(4)*b(2)-d(2)*b(4))+2*zeta(4)*(d(2)*b(1)-d(1)*b(2))+2*zeta(2)*(d(1)*b(4)-d(4)*b(1));

            wfz(1)=2*zeta(3)*(b(1)*c(2)-b(2)*c(1))+2*zeta(1)*(b(2)*c(3)-b(3)*c(2))+2*zeta(2)*(b(3)*c(1)-b(1)*c(3));
            wfz(2)=2*zeta(3)*(b(2)*c(4)-b(4)*c(2))+2*zeta(2)*(b(4)*c(3)-b(3)*c(4))+2*zeta(4)*(b(3)*c(2)-b(2)*c(3));
            wfz(3)=2*zeta(1)*(b(3)*c(4)-b(4)*c(3))+2*zeta(2)*(b(4)*c(1)-b(1)*c(4))+2*zeta(4)*(b(1)*c(3)-b(3)*c(1));
            wfz(4)=2*zeta(1)*(b(4)*c(2)-b(2)*c(4))+2*zeta(4)*(b(2)*c(1)-b(1)*c(2))+2*zeta(2)*(b(1)*c(4)-b(4)*c(1));

            wx(1)=zeta(1)*b(2)-zeta(2)*b(1);      wy(1)=zeta(1)*c(2)-zeta(2)*c(1);      wz(1)=zeta(1)*d(2)-zeta(2)*d(1); 
            wx(2)=zeta(1)*b(3)-zeta(3)*b(1);      wy(2)=zeta(1)*c(3)-zeta(3)*c(1);      wz(2)=zeta(1)*d(3)-zeta(3)*d(1);
            wx(3)=zeta(1)*b(4)-zeta(4)*b(1);      wy(3)=zeta(1)*c(4)-zeta(4)*c(1);      wz(3)=zeta(1)*d(4)-zeta(4)*d(1);
            wx(4)=zeta(2)*b(3)-zeta(3)*b(2);      wy(4)=zeta(2)*c(3)-zeta(3)*c(2);      wz(4)=zeta(2)*d(3)-zeta(3)*d(2);
            wx(5)=zeta(2)*b(4)-zeta(4)*b(2);      wy(5)=zeta(2)*c(4)-zeta(4)*c(2);      wz(5)=zeta(2)*d(4)-zeta(4)*d(2);
            wx(6)=zeta(3)*b(4)-zeta(4)*b(3);      wy(6)=zeta(3)*c(4)-zeta(4)*c(3);      wz(6)=zeta(3)*d(4)-zeta(4)*d(3);
             
            wx=wx.*ll;wy=wy.*ll;wz=wz.*ll;wfx=wfx.*ss;wfy=wfy.*ss;wfz=wfz.*ss;

            xp=xx(ii,jj);yp=yy(ii,jj);zp=zz;dfi=zeros(4,1);dei=zeros(6,1);
            for mm=1:6,if(edges(mm).UknownIndex==0),dei(mm)=1;elseif(sign(edges(mm).UknownIndex)<0),dei(mm)=-1;else,dei(mm)=1;end,end
            for mm=1:4,if(facets(mm).UknownIndex==0),dfi(mm)=1;elseif(sign(facets(mm).UknownIndex)<0),dfi(mm)=-1;else,dfi(mm)=1;end,end
            switch Axis
                case "x"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*xp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*xp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*xp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);
                                                Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*xp);
                                                Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*xp);
                                                Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*xp);
                    end
                case "y"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*yp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*yp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*yp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*yp);
                                                         Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*yp);
                                                         Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*yp);
                    end
                case "z"
                    for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*sgn*dfi(kk)*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*zp);
                                                          By(ii,jj)=By(ii,jj)+wfy(kk)*sgn*dfi(kk)*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*zp);
                                                          Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*sgn*dfi(kk)*FacetVector(facets(kk).Index)*exp(-1i*EigenValue*zp);
                    end
                    for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*zp);
                                                         Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*zp);
                                                         Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index)*exp(-1i*EigenValue*zp);
                    end
            end
         end
    end
    rEx=real(Ex);[ii,jj]=find(abs(rEx)<1000*eps);rEx(ii,jj)=0;rEy=real(Ey);[ii,jj]=find(abs(rEy)<1000*eps);rEy(ii,jj)=0;rEz=real(Ez);[ii,jj]=find(abs(rEz)<1000*eps);rEz(ii,jj)=0;
    iEx=imag(Ex);[ii,jj]=find(abs(iEx)<1000*eps);iEx(ii,jj)=0;iEy=imag(Ey);[ii,jj]=find(abs(iEy)<1000*eps);iEy(ii,jj)=0;iEz=imag(Ez);[ii,jj]=find(abs(iEz)<1000*eps);iEz(ii,jj)=0;
    rBx=real(Bx);[ii,jj]=find(abs(rBx)<1000*eps);rBx(ii,jj)=0;rBy=real(By);[ii,jj]=find(abs(rBy)<1000*eps);rBy(ii,jj)=0;rBz=real(Bz);[ii,jj]=find(abs(rBz)<1000*eps);rBz(ii,jj)=0;
    iBx=imag(Bx);[ii,jj]=find(abs(iBx)<1000*eps);iBx(ii,jj)=0;iBy=imag(By);[ii,jj]=find(abs(iBy)<1000*eps);iBy(ii,jj)=0;iBz=imag(Bz);[ii,jj]=find(abs(iBz)<1000*eps);iBz(ii,jj)=0;
    [cmap]=buildcmap('cbkry');figure;
    subplot(1,3,1);pcolor(xx,yy,rEx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(xx,yy,rEy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(xx,yy,rEz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Electric Field Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(xx,yy,rBx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(xx,yy,rBy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(xx,yy,rBz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Magnetic Flux Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(xx,yy,iEx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(xx,yy,iEy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(xx,yy,iEz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    figure;
    subplot(1,3,1);pcolor(xx,yy,iBx);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(xx,yy,iBy);axis tight;axis equal;shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(xx,yy,iBz);axis tight;axis equal;shading interp;colormap(cmap);title('z Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
end
