function [] = Plot_EH(TModel,BoundaryIndices,X,Nh,Nv)
    if(isscalar(BoundaryIndices)),boundary=TModel.Boundaries(BoundaryIndices);
                                  vertices=TModel.Vertices(boundary.Vertices);
                                  Xs=[vertices.X];Ys=[vertices.Y];Zs=[vertices.Z];
                                  Xmax=max(Xs);Xmin=min(Xs);Ymax=max(Ys);Ymin=min(Ys);Zmax=max(Zs);Zmin=min(Zs);
    else,vertices=[];
        for ii=1:numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(ii));
            temp=vertices;vertices=zeros(numel(boundary.Vertices)+numel(temp),1);vertices(1:numel(temp))=temp;vertices(numel(temp)+1:end)=boundary.Vertices;
        end
        vertices=TModel.Vertices(vertices);
        Xs=[vertices.X];Ys=[vertices.Y];Zs=[vertices.Z];
        Xmax=max(Xs);Xmin=min(Xs);Ymax=max(Ys);Ymin=min(Ys);Zmax=max(Zs);Zmin=min(Zs);
        boundary=TModel.Boundaries(BoundaryIndices(1));Nh=100;Nv=100;
    end
    switch abs(boundary.Axis)
        case 1,yy=linspace(Ymin,Ymax,Nh);zz=linspace(Zmin,Zmax,Nv);Plot_EH_x(TModel,BoundaryIndices,X,yy,zz);
        case 2,xx=linspace(Xmin,Xmax,Nh);zz=linspace(Zmin,Zmax,Nv);Plot_EH_y(TModel,BoundaryIndices,X,xx,zz);
        case 3,xx=linspace(Xmin,Xmax,Nh);yy=linspace(Ymin,Ymax,Nv);Plot_EH_z(TModel,BoundaryIndices,X,xx,yy);
    end
end
%==========================================================================
function [] = Plot_EH_x(TModel,BoundaryIndices,X,yy,zz),Ny=numel(yy);Nz=numel(zz);Ex=zeros(Ny,Nz);Ey=zeros(Ny,Nz);Ez=zeros(Ny,Nz);Hx=zeros(Ny,Nz);Hy=zeros(Ny,Nz);Hz=zeros(Ny,Nz);
   for ib = numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(ib));elements=TModel.Facets(boundary.Facets);BC=[elements.Barycenter];ybc=BC(2,:);zbc=BC(3,:);
       for ii=1:numel(yy)
           for jj=1:numel(zz),ic=0;ie=0;dist=(yy(ii)-ybc).^2 +(zz(jj)-zbc).^2;[~,Is]=sort(dist);
               while ie==0,ic=ic+1;element=elements(Is(ic));
                   edges=[TModel.Edges(element.Edges)];vertices=[TModel.Vertices(element.Vertices)];y=[vertices.Y];z=[vertices.Z];
                   De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');
                   a(1)=(y(1)*z(2)-y(2)*z(1))/De;a(2)=(y(2)*z(3)-y(3)*z(2))/De;a(3)=(y(3)*z(1)-y(1)*z(3))/De;
                   b(1)=(z(1)-z(2))/De;b(2)=(z(2)-z(3))/De;b(3)=(z(3)-z(1))/De;
                   c(1)=(y(2)-y(1))/De;c(2)=(y(3)-y(2))/De;c(3)=(y(1)-y(3))/De;
                   zeta(1)=a(1)+b(1)*yy(ii)+c(1)*zz(jj);
                   zeta(2)=a(2)+b(2)*yy(ii)+c(2)*zz(jj);
                   zeta(3)=a(3)+b(3)*yy(ii)+c(3)*zz(jj);
                if (zeta(1)>0 || abs(zeta(1))<1000*eps) && (zeta(2)>0 || abs(zeta(2))<1000*eps) && (zeta(3)>0 || abs(zeta(3))<1000*eps),ie=Is(ic);end
                if(ic==numel(Is)),ie=Is(1);end
               end,element=elements(ie);
               edges=[TModel.Edges(element.Edges)];vertices=[TModel.Vertices(element.Vertices)];y=[vertices.Y];z=[vertices.Z];
               De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');
               a(1)=(y(1)*z(2)-y(2)*z(1))/De;a(2)=(y(2)*z(3)-y(3)*z(2))/De;a(3)=(y(3)*z(1)-y(1)*z(3))/De;
               b(1)=(z(1)-z(2))/De;b(2)=(z(2)-z(3))/De;b(3)=(z(3)-z(1))/De;
               c(1)=(y(2)-y(1))/De;c(2)=(y(3)-y(2))/De;c(3)=(y(1)-y(3))/De;
               zeta(1)=a(1)+b(1)*yy(ii)+c(1)*zz(jj);
               zeta(2)=a(2)+b(2)*yy(ii)+c(2)*zz(jj);
               zeta(3)=a(3)+b(3)*yy(ii)+c(3)*zz(jj);
               wy(1)=zeta(1)*b(2)-zeta(2)*b(1);wz(1)=zeta(1)*c(2)-zeta(2)*c(1);
               wy(2)=zeta(2)*b(3)-zeta(3)*b(2);wz(2)=zeta(2)*c(3)-zeta(3)*c(2);
               wy(3)=zeta(3)*b(1)-zeta(1)*b(3);wz(3)=zeta(3)*c(1)-zeta(1)*c(3);
               for nn=1:3,edge=edges(nn);vertex=vertices(nn);
                      if(edge.IndexE~=0),si=element.EdgeSigns(nn);Ey(ii,jj)=Ey(ii,jj)+wy(nn)*si*X(edge.IndexE);
                                                                  Ez(ii,jj)=Ez(ii,jj)+wz(nn)*si*X(edge.IndexE);
                      end
                      if(edge.IndexH~=0),si=element.EdgeSigns(nn);Hy(ii,jj)=Hy(ii,jj)+wy(nn)*si*X(edge.IndexH);
                                                                  Hz(ii,jj)=Hz(ii,jj)+wz(nn)*si*X(edge.IndexH);
                      end
                      if(vertex.IndexE~=0),Ex(ii,jj)=Ex(ii,jj)+zeta(nn)*X(vertex.IndexE);end
                      if(vertex.IndexH~=0),Hx(ii,jj)=Hx(ii,jj)+zeta(nn)*X(vertex.IndexH);end
               end
           end
       end
   end
    [cmap]=buildcmap('cbkry');figure;title("Real Fields");
    subplot(2,3,1);pcolor(yy,zz,real(Ex));shading interp;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(yy,zz,real(Ey));shading interp;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(yy,zz,real(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(yy,zz,real(Hx));shading interp;colormap(cmap);title('x Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(yy,zz,real(Hy));shading interp;colormap(cmap);title('y Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(yy,zz,real(Hz));shading interp;colormap(cmap);title('z Component Magnetic Field');colorbar;axis tight;axis equal;
    
    [cmap]=buildcmap('cbkry');figure;title("Imaginary Fields");
    subplot(2,3,1);pcolor(yy,zz,imag(Ex));shading interp;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(yy,zz,imag(Ey));shading interp;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(yy,zz,imag(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(yy,zz,imag(Hx));shading interp;colormap(cmap);title('x Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(yy,zz,imag(Hy));shading interp;colormap(cmap);title('y Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(yy,zz,imag(Hz));shading interp;colormap(cmap);title('z Component Magnetic Field');colorbar;axis tight;axis equal;
end
function [] = Plot_EH_y(TModel,BoundaryIndices,X,xx,zz),Ny=numel(xx);Nz=numel(zz);Ex=zeros(Ny,Nz);Ey=zeros(Ny,Nz);Ez=zeros(Ny,Nz);Hx=zeros(Ny,Nz);Hy=zeros(Ny,Nz);Hz=zeros(Ny,Nz);
   for ib = numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(ib));elements=TModel.Facets(boundary.Facets);BC=[elements.Barycenter];xbc=BC(1,:);zbc=BC(3,:);
       for ii=1:numel(xx)
           for jj=1:numel(zz),ic=0;ie=0;dist=(xx(ii)-xbc).^2 +(zz(jj)-zbc).^2;[~,Is]=sort(dist);
               while ie==0,ic=ic+1;element=elements(Is(ic));
                   vertices=[TModel.Vertices(element.Vertices)];x=[vertices.X];y=[vertices.Y];z=[vertices.Z];
                   De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');
                   a(1)=(x(1)*z(2)-x(2)*z(1))/De;a(2)=(x(2)*z(3)-x(3)*z(2))/De;a(3)=(x(3)*z(1)-x(1)*z(3))/De;
                   b(1)=(z(1)-z(2))/De;b(2)=(z(2)-z(3))/De;b(3)=(z(3)-z(1))/De;
                   c(1)=(x(2)-x(1))/De;c(2)=(x(3)-x(2))/De;c(3)=(x(1)-x(3))/De;
                   zeta(1)=a(1)+b(1)*xx(ii)+c(1)*zz(jj);
                   zeta(2)=a(2)+b(2)*xx(ii)+c(2)*zz(jj);
                   zeta(3)=a(3)+b(3)*xx(ii)+c(3)*zz(jj);
                if (zeta(1)>0 || abs(zeta(1))<1000*eps) && (zeta(2)>0 || abs(zeta(2))<1000*eps) && (zeta(3)>0 || abs(zeta(3))<1000*eps),ie=Is(ic);end
                if(ic==numel(Is)),ie=Is(1);end
               end,element=elements(ie);
               edges=[TModel.Edges(element.Edges)];vertices=[TModel.Vertices(element.Vertices)];x=[vertices.X];y=[vertices.Y];z=[vertices.Z];
               De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');
               a(1)=(x(1)*z(2)-x(2)*z(1))/De;a(2)=(x(2)*z(3)-x(3)*z(2))/De;a(3)=(x(3)*z(1)-x(1)*z(3))/De;
               b(1)=(z(1)-z(2))/De;b(2)=(z(2)-z(3))/De;b(3)=(z(3)-z(1))/De;
               c(1)=(x(2)-x(1))/De;c(2)=(x(3)-x(2))/De;c(3)=(x(1)-x(3))/De;
               zeta(1)=a(1)+b(1)*xx(ii)+c(1)*zz(jj);
               zeta(2)=a(2)+b(2)*xx(ii)+c(2)*zz(jj);
               zeta(3)=a(3)+b(3)*xx(ii)+c(3)*zz(jj);
               wx(1)=zeta(1)*b(2)-zeta(2)*b(1);wz(1)=zeta(1)*c(2)-zeta(2)*c(1);
               wx(2)=zeta(2)*b(3)-zeta(3)*b(2);wz(2)=zeta(2)*c(3)-zeta(3)*c(2);
               wx(3)=zeta(3)*b(1)-zeta(1)*b(3);wz(3)=zeta(3)*c(1)-zeta(1)*c(3);
               for nn=1:3,edge=edges(nn);vertex=vertices(nn);
                      if(edge.IndexE~=0),si=element.EdgeSigns(nn);Ex(ii,jj)=Ex(ii,jj)+wx(nn)*si*X(edge.IndexE);
                                                                  Ez(ii,jj)=Ez(ii,jj)+wz(nn)*si*X(edge.IndexE);
                      end
                      if(edge.IndexH~=0),si=element.EdgeSigns(nn);Hx(ii,jj)=Hx(ii,jj)+wy(nn)*si*X(edge.IndexH);
                                                                  Hz(ii,jj)=Hz(ii,jj)+wz(nn)*si*X(edge.IndexH);
                      end
                      if(vertex.IndexE~=0),Ey(ii,jj)=Ey(ii,jj)+zeta(nn)*X(vertex.IndexE);end
                      if(vertex.IndexH~=0),Hy(ii,jj)=Hy(ii,jj)+zeta(nn)*X(vertex.IndexH);end
               end
           end
       end
   end
    [cmap]=buildcmap('cbkry');figure;title("Real Fields");
    subplot(2,3,1);pcolor(xx,zz,real(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(xx,zz,real(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(xx,zz,real(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(xx,zz,real(Hx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(xx,zz,real(Hy));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(xx,zz,real(Hz));shading interp;colormap(cmap);title('z Component Magnetic Field');colorbar;axis tight;axis equal;
    
    [cmap]=buildcmap('cbkry');figure;title("Imaginary Fields");
    subplot(2,3,1);pcolor(xx,zz,imag(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(xx,zz,imag(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(xx,zz,imag(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(xx,zz,imag(Hx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(xx,zz,imag(Hy));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(xx,zz,imag(Hz));shading interp;colormap(cmap);title('z Component Magnetic Field');colorbar;axis tight;axis equal;
end
function [] = Plot_EH_z(TModel,BoundaryIndices,X,xx,yy),Nx=numel(xx);Ny=numel(yy);Ex=zeros(Nx,Ny);Ey=zeros(Nx,Ny);Ez=zeros(Nx,Ny);Hx=zeros(Nx,Ny);Hy=zeros(Nx,Ny);Hz=zeros(Nx,Ny);
   for ib = numel(BoundaryIndices),boundary=TModel.Boundaries(BoundaryIndices(ib));elements=TModel.Facets(boundary.Facets);BC=[elements.Barycenter];xbc=BC(1,:);ybc=BC(2,:);
       for ii=1:numel(xx)
           for jj=1:numel(yy),ic=0;ie=0;dist=(xx(ii)-xbc).^2 +(yy(jj)-ybc).^2;[~,Is]=sort(dist);
               while ie==0,ic=ic+1;element=elements(Is(ic));
                   vertices=[TModel.Vertices(element.Vertices)];x=[vertices.X];y=[vertices.Y];z=[vertices.Z];
                   De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');
                   a(1)=(x(1)*y(2)-x(2)*y(1))/De;a(2)=(x(2)*y(3)-x(3)*y(2))/De;a(3)=(x(3)*y(1)-x(1)*y(3))/De;
                   b(1)=(y(1)-y(2))/De;b(2)=(y(2)-y(3))/De;b(3)=(y(3)-y(1))/De;
                   c(1)=(x(2)-x(1))/De;c(2)=(x(3)-x(2))/De;c(3)=(x(1)-x(3))/De;
                   zeta(1)=a(1)+b(1)*xx(ii)+c(1)*yy(jj);
                   zeta(2)=a(2)+b(2)*xx(ii)+c(2)*yy(jj);
                   zeta(3)=a(3)+b(3)*xx(ii)+c(3)*yy(jj);
                if (zeta(1)>0 || abs(zeta(1))<1000*eps) && (zeta(2)>0 || abs(zeta(2))<1000*eps) && (zeta(3)>0 || abs(zeta(3))<1000*eps),ie=Is(ic);end
                if(ic==numel(Is)),ie=Is(1);end
               end,element=elements(ie);
               edges=[TModel.Edges(element.Edges)];vertices=[TModel.Vertices(element.Vertices)];x=[vertices.X];y=[vertices.Y];z=[vertices.Z];
               De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');
               a(1)=(x(1)*y(2)-x(2)*y(1))/De;a(2)=(x(2)*y(3)-x(3)*y(2))/De;a(3)=(x(3)*y(1)-x(1)*y(3))/De;
               b(1)=(y(1)-y(2))/De;b(2)=(y(2)-y(3))/De;b(3)=(y(3)-y(1))/De;
               c(1)=(x(2)-x(1))/De;c(2)=(x(3)-x(2))/De;c(3)=(x(1)-x(3))/De;
               zeta(1)=a(1)+b(1)*xx(ii)+c(1)*yy(jj);
               zeta(2)=a(2)+b(2)*xx(ii)+c(2)*yy(jj);
               zeta(3)=a(3)+b(3)*xx(ii)+c(3)*yy(jj);
               wx(1)=zeta(1)*b(2)-zeta(2)*b(1);wy(1)=zeta(1)*c(2)-zeta(2)*c(1);
               wx(2)=zeta(2)*b(3)-zeta(3)*b(2);wy(2)=zeta(2)*c(3)-zeta(3)*c(2);
               wx(3)=zeta(3)*b(1)-zeta(1)*b(3);wy(3)=zeta(3)*c(1)-zeta(1)*c(3);
               for nn=1:3,edge=edges(nn);vertex=vertices(nn);
                      if(edge.IndexE~=0),si=element.EdgeSigns(nn);Ex(ii,jj)=Ex(ii,jj)+wx(nn)*si*X(edge.IndexE);
                                                                  Ey(ii,jj)=Ey(ii,jj)+wy(nn)*si*X(edge.IndexE);
                      end
                      if(edge.IndexH~=0),si=element.EdgeSigns(nn);Hx(ii,jj)=Hx(ii,jj)+wx(nn)*si*X(edge.IndexH);
                                                                  Hy(ii,jj)=Hy(ii,jj)+wy(nn)*si*X(edge.IndexH);
                      end
                      if(vertex.IndexE~=0),Ez(ii,jj)=Ez(ii,jj)+zeta(nn)*X(vertex.IndexE);end
                      if(vertex.IndexH~=0),Hz(ii,jj)=Hz(ii,jj)+zeta(nn)*X(vertex.IndexH);end
               end
           end
       end
   end
    [cmap]=buildcmap('cbkry');figure;title("Real Fields");
    subplot(2,3,1);pcolor(xx,yy,real(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(xx,yy,real(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(xx,yy,real(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(xx,yy,real(Hx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(xx,yy,real(Hy));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(xx,yy,real(Hz));shading interp;colormap(cmap);title('z Component Magnetic Field');colorbar;axis tight;axis equal;
    
    [cmap]=buildcmap('cbkry');figure;title("Imaginary Fields");
    subplot(2,3,1);pcolor(xx,yy,imag(Ex));shading interp;hold on;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,2);pcolor(xx,yy,imag(Ey));shading interp;hold on;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,3);pcolor(xx,yy,imag(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
    subplot(2,3,4);pcolor(xx,yy,imag(Hx));shading interp;hold on;colormap(cmap);title('x Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,5);pcolor(xx,yy,imag(Hy));shading interp;hold on;colormap(cmap);title('y Component Magnetic Field');colorbar;axis tight;axis equal;
    subplot(2,3,6);pcolor(xx,yy,imag(Hz));shading interp;colormap(cmap);title('z Component Magnetic Field');colorbar;axis tight;axis equal;
end