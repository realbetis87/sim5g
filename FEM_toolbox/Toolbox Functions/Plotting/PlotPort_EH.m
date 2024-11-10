function [] = PlotPort_EH(TModel,BoundaryIndices,X)
    if(isscalar(BoundaryIndices)),boundary=TModel.Boundaries(BoundaryIndices);else,boundary=TModel.Boundaries(BoundaryIndices(1));end
    switch abs(boundary.Axis)
        case 1,PlotEH_x(TModel,BoundaryIndices,X);
        case 2,PlotEH_y(TModel,BoundaryIndices,X);
        case 3,PlotEH_z(TModel,BoundaryIndices,X);
    end
end
function [] = PlotEH_x2(TModel,BoundaryIndices,X)
    Nh=100;Nv=50;Ex=zeros(Nh,Nv);Ey=zeros(Nh,Nv);Ez=zeros(Nh,Nv);Hx=zeros(Nh,Nv);Hy=zeros(Nh,Nv);Hz=zeros(Nh,Nv);
    boundary=TModel.Boundaries(BoundaryIndices);
    vertices=TModel.Vertices(boundary.Vertices);elements=TModel.Facets(boundary.Facets);
    ymin=min([vertices.Y]);ymax=max([vertices.Y]);zmin=min([vertices.Z]);zmax=max([vertices.Z]);
    yy=linspace(ymin,ymax,Nh);zz=linspace(zmin,zmax,Nv);
    for ii=1:Nh
        for jj=1:Nv
            for kk=1:numel(boundary.Facets),element=TModel.Facets(boundary.Facets(kk));
                edges=[TModel.Edges(element.Edges)];vertices=[TModel.Vertices(element.Vertices)];x=[vertices.X];y=[vertices.Y];z=[vertices.Z];
                De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');
                a(1)=(y(1)*z(2)-y(2)*z(1))/De;a(2)=(y(2)*z(3)-y(3)*z(2))/De;a(3)=(y(3)*z(1)-y(1)*z(3))/De;
                b(1)=(z(1)-z(2))/De;b(2)=(z(2)-z(3))/De;b(3)=(z(3)-z(1))/De;
                c(1)=(y(2)-y(1))/De;c(2)=(y(3)-y(2))/De;c(3)=(y(1)-y(3))/De;
                zeta(1)=a(1)+b(1)*yy(ii)+c(1)*zz(jj);
                zeta(2)=a(2)+b(2)*yy(ii)+c(2)*zz(jj);
                zeta(3)=a(3)+b(3)*yy(ii)+c(3)*zz(jj);
                if (zeta(1)>0 || abs(zeta(1))<1000*eps) && (zeta(2)>0 || abs(zeta(2))<1000*eps) && (zeta(3)>0 || abs(zeta(3))<1000*eps)
                  wy(1)=zeta(1)*b(2)-zeta(2)*b(1);wz(1)=zeta(1)*c(2)-zeta(2)*c(1);
                  wy(2)=zeta(2)*b(3)-zeta(3)*b(2);wz(2)=zeta(2)*c(3)-zeta(3)*c(2);
                  wy(3)=zeta(3)*b(1)-zeta(1)*b(3);wz(3)=zeta(3)*c(1)-zeta(1)*c(3);
                   for mm=1:3,edge=edges(mm);vertex=vertices(mm);
                      if(edge.IndexE~=0),si=element.EdgeSigns(mm);Ey(ii,jj)=Ey(ii,jj)+wy(mm)*si*X(edge.IndexE);Ez(ii,jj)=Ez(ii,jj)+wz(mm)*si*X(edge.IndexE);end
                      if(edge.IndexH~=0),si=element.EdgeSigns(mm);Hy(ii,jj)=Hy(ii,jj)+wy(mm)*si*X(edge.IndexH);Hz(ii,jj)=Hz(ii,jj)+wz(mm)*si*X(edge.IndexH);end
                      if(vertex.IndexE~=0),Ex(ii,jj)=Ex(ii,jj)+zeta(mm)*X(vertex.IndexE);end
                      if(vertex.IndexH~=0),Hx(ii,jj)=Hx(ii,jj)+zeta(mm)*X(vertex.IndexH);end
                  end
                end
            end
        end
    end
     [cmap]=buildcmap('cbkry');figure;title("Real Fields");
    subplot(2,3,1);pcolor(yy,zz,real(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(yy,zz,real(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(yy,zz,real(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(yy,zz,real(Hx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(yy,zz,real(Hy));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(yy,zz,real(Hz));shading interp;colormap(cmap);title('z Component Magnetic Field');colorbar;axis tight;axis equal;
    
    [cmap]=buildcmap('cbkry');figure;title("Imaginary Fields");
    subplot(2,3,1);pcolor(yy,zz,imag(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(yy,zz,imag(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(yy,zz,imag(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(yy,zz,imag(Hx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(yy,zz,imag(Hy));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(yy,zz,imag(Hz));shading interp;colormap(cmap);title('z Component Magnetic Field');colorbar;axis tight;axis equal;
end
function [] = PlotEH_x(TModel,BoundaryIndices,X)
    NE=0;Nv=zeros(numel(BoundaryIndices),1);
    for kk=1:numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(kk));NE=NE+numel(boundary.Facets);
        verts=[TModel.Vertices(boundary.Vertices)];Nv(kk)=max([verts.Index2D]);
    end,Trian=zeros(NE,3);NV=max(Nv);
    Ex=zeros(NE,1);Ey=zeros(NE,1);Ez=zeros(NE,1);Hx=zeros(NE,1);Hy=zeros(NE,1);Hz=zeros(NE,1);Ys=zeros(NV,1);Zs=zeros(NV,1);counter=0;
    YY=zeros(NE,1);ZZ=zeros(NE,1);
    for kk=1:numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=TModel.Facets(boundary.Facets(ie));edges=[TModel.Edges(element.Edges)];vertices=[TModel.Vertices(element.Vertices)];x=[vertices.X];y=[vertices.Y];z=[vertices.Z];
                  counter=counter+1;
                  Ys(vertices(1).Index2D)=y(1);Ys(vertices(2).Index2D)=y(2);Ys(vertices(3).Index2D)=y(3);
                  Zs(vertices(1).Index2D)=z(1);Zs(vertices(2).Index2D)=z(2);Zs(vertices(3).Index2D)=z(3);
                  Trian(counter,1)=vertices(1).Index2D;Trian(counter,2)=vertices(2).Index2D;Trian(counter,3)=vertices(3).Index2D;
                  YY(counter)=element.Barycenter(2);ZZ(counter)=element.Barycenter(3);
                  De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');
                  a(1)=(y(1)*z(2)-y(2)*z(1))/De;a(2)=(y(2)*z(3)-y(3)*z(2))/De;a(3)=(y(3)*z(1)-y(1)*z(3))/De;
                  b(1)=(z(1)-z(2))/De;b(2)=(z(2)-z(3))/De;b(3)=(z(3)-z(1))/De;
                  c(1)=(y(2)-y(1))/De;c(2)=(y(3)-y(2))/De;c(3)=(y(1)-y(3))/De;
                  zeta(1)=a(1)+b(1)*element.Barycenter(2)+c(1)*element.Barycenter(3);
                  zeta(2)=a(2)+b(2)*element.Barycenter(2)+c(2)*element.Barycenter(3);
                  zeta(3)=a(3)+b(3)*element.Barycenter(2)+c(3)*element.Barycenter(3);
                  wy(1)=zeta(1)*b(2)-zeta(2)*b(1);wz(1)=zeta(1)*c(2)-zeta(2)*c(1);
                  wy(2)=zeta(2)*b(3)-zeta(3)*b(2);wz(2)=zeta(2)*c(3)-zeta(3)*c(2);
                  wy(3)=zeta(3)*b(1)-zeta(1)*b(3);wz(3)=zeta(3)*c(1)-zeta(1)*c(3);
                  for ii=1:3,edge=edges(ii);vertex=vertices(ii);
                      if(edge.IndexE~=0),si=element.EdgeSigns(ii);Ey(counter)=Ey(counter)+wy(ii)*si*X(edge.IndexE);
                                                                  Ez(counter)=Ez(counter)+wz(ii)*si*X(edge.IndexE);
                      end
                      if(edge.IndexH~=0),si=element.EdgeSigns(ii);Hy(counter)=Hy(counter)+wy(ii)*si*X(edge.IndexH);
                                                                  Hz(counter)=Hz(counter)+wz(ii)*si*X(edge.IndexH);
                      end
                      if(vertex.IndexE~=0),Ex(counter)=Ex(counter)+zeta(ii)*X(vertex.IndexE);end
                      if(vertex.IndexH~=0),Hx(counter)=Hx(counter)+zeta(ii)*X(vertex.IndexH);end
                  end

            end
    end
    [H,V] = meshgrid(YY,ZZ);rEx=griddata(YY,ZZ,real(Ex),H,V,'v4');iEx=griddata(YY,ZZ,imag(Ex),H,V,'v4');
                            rEy=griddata(YY,ZZ,real(Ey),H,V,'v4');iEy=griddata(YY,ZZ,imag(Ey),H,V,'v4');
                            rEz=griddata(YY,ZZ,real(Ez),H,V,'v4');iEz=griddata(YY,ZZ,imag(Ez),H,V,'v4');
                            rHx=griddata(YY,ZZ,real(Hx),H,V,'v4');iHx=griddata(YY,ZZ,imag(Hx),H,V,'v4');
                            rHy=griddata(YY,ZZ,real(Hy),H,V,'v4');iHy=griddata(YY,ZZ,imag(Hy),H,V,'v4');
                            rHz=griddata(YY,ZZ,real(Hz),H,V,'v4');iHz=griddata(YY,ZZ,imag(Hz),H,V,'v4');
    [cmap]=buildcmap('cbkry');figure;
    subplot(2,3,1);pcolor(YY,ZZ,rEx);shading interp;hold on;colormap(cmap);title('x Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(YY,ZZ,rEy);shading interp;hold on;colormap(cmap);title('y Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(YY,ZZ,rEz);shading interp;colormap(cmap);title('z Component Electric Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(YY,ZZ,rHx);shading interp;hold on;colormap(cmap);title('x Component Magnetic Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(YY,ZZ,rHy);shading interp;hold on;colormap(cmap);title('y Component Magnetic Field Real Part');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(YY,ZZ,rHz);shading interp;colormap(cmap);title('z Component Magnetic Field Real Part');colorbar;axis tight;axis equal;
    figure;
    subplot(2,3,1);pcolor(YY,ZZ,iEx);shading interp;hold on;colormap(cmap);title('x Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(YY,ZZ,iEy);shading interp;hold on;colormap(cmap);title('y Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(YY,ZZ,iEz);shading interp;colormap(cmap);title('z Component Electric Field Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(YY,ZZ,iHx);shading interp;hold on;colormap(cmap);title('x Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(YY,ZZ,iHy);shading interp;hold on;colormap(cmap);title('y Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(YY,ZZ,iHz);shading interp;colormap(cmap);title('z Component Magnetic Flux Imaginary Part');colorbar;axis tight;axis equal;


end
function [] = PlotEH_y(TModel,BoundaryIndices,X)
    NE=0;for kk=1:numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(kk));NE=NE+numel(boundary.Facets);end
    Ex=zeros(NE,1);Ey=zeros(NE,1);Ez=zeros(NE,1);Hx=zeros(NE,1);Hy=zeros(NE,1);Hz=zeros(NE,1);Xs=zeros(NE,1);Zs=zeros(NE,1);counter=0;
    for kk=1:numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=TModel.Facets(boundary.Facets(ie));edges=[TModel.Edges(element.Edges)];vertices=[TModel.Vertices(element.Vertices)];x=[vertices.X];y=[vertices.Y];z=[vertices.Z];
                  counter=counter+1;Xs(counter)=element.Barycenter(1);Zs(counter)=element.Barycenter(3);
                  De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');Ae=abs(De)/2;
                  a(1)=(x(1)*z(2)-x(2)*z(1))/De;a(2)=(x(2)*z(3)-x(3)*z(2))/De;a(3)=(x(3)*z(1)-x(1)*z(3))/De;
                  b(1)=(z(1)-z(2))/De;b(2)=(z(2)-z(3))/De;b(3)=(z(3)-z(1))/De;
                  c(1)=(x(2)-x(1))/De;c(2)=(x(3)-x(2))/De;c(3)=(x(1)-x(3))/De;
                  zeta(1)=a(1)+b(1)*element.Barycenter(1)+c(1)*element.Barycenter(3);
                  zeta(2)=a(2)+b(2)*element.Barycenter(1)+c(2)*element.Barycenter(3);
                  zeta(3)=a(3)+b(3)*element.Barycenter(1)+c(3)*element.Barycenter(3);
                  wx(1)=zeta(1)*b(2)-zeta(2)*b(1);wz(1)=zeta(1)*c(2)-zeta(2)*c(1);
                  wx(2)=zeta(2)*b(3)-zeta(3)*b(2);wz(2)=zeta(2)*c(3)-zeta(3)*c(2);
                  wx(3)=zeta(3)*b(1)-zeta(1)*b(3);wz(3)=zeta(3)*c(1)-zeta(1)*c(3);
                  for ii=1:3,edge=edges(ii);vertex=vertices(ii);
                      if(edge.IndexE~=0),si=element.EdgeSigns(ii);Ex(counter)=Ex(counter)+wx(ii)*si*X(edge.IndexE);Ez(counter)=Ez(counter)+wz(ii)*si*X(edge.IndexE);end
                      if(edge.IndexH~=0),si=element.EdgeSigns(ii);Hx(counter)=Hx(counter)+wx(ii)*si*X(edge.IndexE);Hz(counter)=Hz(counter)+wz(ii)*si*X(edge.IndexH);end
                      if(vertex.IndexE~=0),Ey(counter)=Ey(counter)+zeta(ii)*X(vertex.IndexE);end
                      if(vertex.IndexH~=0),Hy(counter)=Hy(counter)+zeta(ii)*X(vertex.IndexH);end
                  end

            end
    end
    figure;
end
function [] = PlotEH_z(TModel,BoundaryIndices,X)
    NE=0;for kk=1:numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(kk));NE=NE+numel(boundary.Facets);end
    Ex=zeros(NE,1);Ey=zeros(NE,1);Ez=zeros(NE,1);Hx=zeros(NE,1);Hy=zeros(NE,1);Hz=zeros(NE,1);Xs=zeros(NE,1);Ys=zeros(NE,1);counter=0;
    for kk=1:numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(kk));
            for ie=1:numel(boundary.Facets),element=TModel.Facets(boundary.Facets(ie));edges=[TModel.Edges(element.Edges)];vertices=[TModel.Vertices(element.Vertices)];x=[vertices.X];y=[vertices.Y];z=[vertices.Z];
                  counter=counter+1;Xs(counter)=element.Barycenter(1);Ys(counter)=element.Barycenter(2);
                  De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');Ae=abs(De)/2;
                  a(1)=(x(1)*y(2)-x(2)*y(1))/De;a(2)=(x(2)*y(3)-x(3)*y(2))/De;a(3)=(x(3)*y(1)-x(1)*y(3))/De;
                  b(1)=(y(1)-y(2))/De;b(2)=(y(2)-y(3))/De;b(3)=(y(3)-y(1))/De;
                  c(1)=(x(2)-x(1))/De;c(2)=(x(3)-x(2))/De;c(3)=(x(1)-x(3))/De;
                  zeta(1)=a(1)+b(1)*element.Barycenter(1)+c(1)*element.Barycenter(2);
                  zeta(2)=a(2)+b(2)*element.Barycenter(1)+c(2)*element.Barycenter(2);
                  zeta(3)=a(3)+b(3)*element.Barycenter(1)+c(3)*element.Barycenter(2);
                  wx(1)=zeta(1)*b(2)-zeta(2)*b(1);wy(1)=zeta(1)*c(2)-zeta(2)*c(1);
                  wx(2)=zeta(2)*b(3)-zeta(3)*b(2);wy(2)=zeta(2)*c(3)-zeta(3)*c(2);
                  wx(3)=zeta(3)*b(1)-zeta(1)*b(3);wy(3)=zeta(3)*c(1)-zeta(1)*c(3);
                  for ii=1:3,edge=edges(ii);vertex=vertices(ii);
                      if(edge.IndexE~=0),si=element.EdgeSigns(ii);Ex(counter)=Ex(counter)+wx(ii)*si*X(edge.IndexE);Ey(counter)=Ey(counter)+wy(ii)*si*X(edge.IndexE);end
                      if(edge.IndexH~=0),si=element.EdgeSigns(ii);Hx(counter)=Hx(counter)+wx(ii)*si*X(edge.IndexE);Hy(counter)=Hy(counter)+wy(ii)*si*X(edge.IndexH);end
                      if(vertex.IndexE~=0),Ez(counter)=Ez(counter)+zeta(ii)*X(vertex.IndexE);end
                      if(vertex.IndexH~=0),Hz(counter)=Hz(counter)+zeta(ii)*X(vertex.IndexH);end
                  end

            end
    end
    figure;
end