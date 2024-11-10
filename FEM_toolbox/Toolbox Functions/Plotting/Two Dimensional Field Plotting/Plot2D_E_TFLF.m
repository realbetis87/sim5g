function [] = Plot2D_E_TFLF(TModel,BoundaryExcitation,Nhor,Nver,X)
    NVertices=numel(BoundaryExcitation.Vertices);NEdges=numel(BoundaryExcitation.Edges);NFacets=numel(BoundaryExcitation.Edges);
    VertexVector_E=zeros(NVertices,1);EdgeVector_E=zeros(NEdges,1);
    for ii=1:NEdges,edge=TModel.Edges(BoundaryExcitation.Edges(ii));if(edge.IndexE~=0),EdgeVector_E(edge.Index2D)=X(edge.IndexE);end,end
    for ii=1:NVertices,vertex=TModel.Vertices(BoundaryExcitation.Vertices(ii));if(vertex.IndexE~=0),VertexVector_E(vertex.Index2D)=X(vertex.IndexE);end,end
    Ex=zeros(Nver,Nhor);Ey=zeros(Nver,Nhor);Ez=zeros(Nver,Nhor);
    boundary=TModel.Boundaries(BoundaryExcitation.BoundaryIndices(1));
    Vertices=[TModel.Vertices(BoundaryExcitation.Vertices)];Xs=[Vertices.X];Ys=[Vertices.Y];Zs=[Vertices.Z];
    Xmax=max(Xs);Xmin=min(Xs);Ymax=max(Ys);Ymin=min(Ys);Zmax=max(Zs);Zmin=min(Zs);
    switch abs(boundary.Axis)
        case 1,y=linspace(Ymin-Ymin*1e-3,Ymax-Ymax*1e-3,Nhor);z=linspace(Zmin-Zmin*1e-3,Zmax-Zmax*1e-3,Nver);[yy,zz]=meshgrid(y,z);Plot2D_E_TFLF_x(TModel,BoundaryExcitation,Ex,Ey,Ez,VertexVector_E,EdgeVector_E,yy,zz);
        case 2,x=linspace(Xmin-Xmin*1e-3,Xmax-Xmax*1e-3,Nhor);z=linspace(Zmin-Zmin*1e-3,Zmax-Zmax*1e-3,Nver);[xx,zz]=meshgrid(x,z);Plot2D_E_TFLF_y(TModel,BoundaryExcitation,Ex,Ey,Ez,VertexVector_E,EdgeVector_E,xx,zz);
        case 3,x=linspace(Xmin-Xmin*1e-3,Xmax-Xmax*1e-3,Nhor);y=linspace(Ymin-Ymin*1e-3,Ymax-Ymax*1e-3,Nver);[xx,yy]=meshgrid(x,y);Plot2D_E_TFLF_z(TModel,BoundaryExcitation,Ex,Ey,Ez,VertexVector_E,EdgeVector_E,xx,yy);
    end
end
%--------------------------------------------------------------------------
function [] = Plot2D_E_TFLF_x(TModel,BoundaryExcitation,Ex,Ey,Ez,VertexVector_E,EdgeVector_E,yy,zz),ci=0;
    elements=[TModel.Facets(BoundaryExcitation.Facets)];BaryCenters=[elements.Barycenter];ybc=BaryCenters(2,:);zbc=BaryCenters(3,:);
    for ii=1:size(yy,1)
        for jj=1:size(zz,2),ic=0;ie=0;dist=(yy(ii,jj)-ybc).^2 +(zz(ii,jj)-zbc).^2;[~,Is]=sort(dist);
            while ie==0,ic=ic+1;element=elements(Is(ic));vertices=[TModel.Vertices(element.Vertices)];y=[vertices.Y];z=[vertices.Z];
                De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');
                bb(1)=(z(2)-z(3))/De;  cc(1)=(y(3)-y(2))/De;  aa(1)=(y(2)*z(3)-y(3)*z(2))/De;
                bb(2)=(z(3)-z(1))/De;  cc(2)=(y(1)-y(3))/De;  aa(2)=(y(3)*z(1)-y(1)*z(3))/De;
                bb(3)=(z(1)-z(2))/De;  cc(3)=(y(2)-y(1))/De;  aa(3)=(y(1)*z(2)-y(2)*z(1))/De;
                %----------------------------------------------------------
                zeta(1)=aa(1)+bb(1)*yy(ii,jj)+cc(1)*zz(ii,jj);
                zeta(2)=aa(2)+bb(2)*yy(ii,jj)+cc(2)*zz(ii,jj);
                zeta(3)=aa(3)+bb(3)*yy(ii,jj)+cc(3)*zz(ii,jj);
                if (zeta(1)>0 || abs(zeta(1))<1000*eps) && (zeta(2)>0 || abs(zeta(2))<1000*eps) && (zeta(3)>0 || abs(zeta(3))<1000*eps),ie=Is(ic);end
                if(ic==numel(Is)),ie=Is(1);ci=ci+1;
                    element=elements(Is(1));vertices=[TModel.Vertices(element.Vertices)];y=[vertices.Y];z=[vertices.Z];
                    De=det([1 y(1) z(1);1 y(2) z(2);1 y(3) z(3);]');
                    bb(1)=(z(2)-z(3))/De;  cc(1)=(y(3)-y(2))/De;  aa(1)=(y(2)*z(3)-y(3)*z(2))/De;
                    bb(2)=(z(3)-z(1))/De;  cc(2)=(y(1)-y(3))/De;  aa(2)=(y(3)*z(1)-y(1)*z(3))/De;
                    bb(3)=(z(1)-z(2))/De;  cc(3)=(y(2)-y(1))/De;  aa(3)=(y(1)*z(2)-y(2)*z(1))/De;
                    %----------------------------------------------------------
                    zeta(1)=aa(1)+bb(1)*yy(ii,jj)+cc(1)*zz(ii,jj);
                    zeta(2)=aa(2)+bb(2)*yy(ii,jj)+cc(2)*zz(ii,jj);
                    zeta(3)=aa(3)+bb(3)*yy(ii,jj)+cc(3)*zz(ii,jj);
                end
            end,edges=[TModel.Edges(element.Edges)];ll=[edges.Length];
            wwy(1)=ll(1)*(zeta(1)*bb(2)-zeta(2)*bb(1));wwz(1)=ll(1)*(zeta(1)*cc(2)-zeta(2)*cc(1));wwy(2)=ll(2)*(zeta(2)*bb(3)-zeta(3)*bb(2));
            wwz(2)=ll(2)*(zeta(2)*cc(3)-zeta(3)*cc(2));wwy(3)=ll(3)*(zeta(3)*bb(1)-zeta(1)*bb(3));wwz(3)=ll(3)*(zeta(3)*cc(1)-zeta(1)*cc(3));
            for kk=1:3
                Ey(ii,jj)=Ey(ii,jj)+element.EdgeSigns(kk)*wwy(kk)*EdgeVector_E(edges(kk).Index2D);
                Ez(ii,jj)=Ez(ii,jj)+element.EdgeSigns(kk)*wwz(kk)*EdgeVector_E(edges(kk).Index2D);
                Ex(ii,jj)=Ex(ii,jj)+zeta(kk)*VertexVector_E(vertices(kk).Index2D);
            end
        end
    end
    [cmap]=buildcmap('cbkry');figure;title("Real Fields");
    subplot(1,3,1);pcolor(yy,zz,real(Ex));shading interp;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(yy,zz,real(Ey));shading interp;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(yy,zz,real(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
    [cmap]=buildcmap('cbkry');figure;title("Imaginary Fields");
    subplot(1,3,1);pcolor(yy,zz,imag(Ex));shading interp;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(yy,zz,imag(Ey));shading interp;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(yy,zz,imag(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
   [cmap]=buildcmap('cbkry');figure;title("Absolute Fields");
    subplot(1,3,1);pcolor(yy,zz,abs(Ex));shading interp;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(yy,zz,abs(Ey));shading interp;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(yy,zz,abs(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
end
function [] = Plot2D_E_TFLF_y(TModel,BoundaryExcitation,Ex,Ey,Ez,VertexVector_E,EdgeVector_E,xx,zz),ci=0;
    elements=[TModel.Facets(BoundaryExcitation.Facets)];BaryCenters=[elements.Barycenter];xbc=BaryCenters(1,:);zbc=BaryCenters(3,:);
    for ii=1:size(xx,1)
        for jj=1:size(zz,2),ic=0;ie=0;dist=(xx(ii,jj)-xbc).^2 +(zz(ii,jj)-zbc).^2;[~,Is]=sort(dist);
            while ie==0,ic=ic+1;element=elements(Is(ic));vertices=[TModel.Vertices(element.Vertices)];x=[vertices.X];z=[vertices.Z];
                De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');
                bb(1)=(z(2)-z(3))/De;  cc(1)=(x(3)-x(2))/De;  aa(1)=(x(2)*z(3)-x(3)*z(2))/De;
                bb(2)=(z(3)-z(1))/De;  cc(2)=(x(1)-x(3))/De;  aa(2)=(x(3)*z(1)-x(1)*z(3))/De;
                bb(3)=(z(1)-z(2))/De;  cc(3)=(x(2)-x(1))/De;  aa(3)=(x(1)*z(2)-x(2)*z(1))/De;
                %----------------------------------------------------------
                zeta(1)=aa(1)+bb(1)*xx(ii,jj)+cc(1)*zz(ii,jj);
                zeta(2)=aa(2)+bb(2)*xx(ii,jj)+cc(2)*zz(ii,jj);
                zeta(3)=aa(3)+bb(3)*xx(ii,jj)+cc(3)*zz(ii,jj);
                if (zeta(1)>0 || abs(zeta(1))<1000*eps) && (zeta(2)>0 || abs(zeta(2))<1000*eps) && (zeta(3)>0 || abs(zeta(3))<1000*eps),ie=Is(ic);end
                if(ic==numel(Is)),ie=Is(1);ci=ci+1;%disp(ci);
                    element=elements(Is(1));vertices=[TModel.Vertices(element.Vertices)];x=[vertices.X];z=[vertices.Z];
                    De=det([1 x(1) z(1);1 x(2) z(2);1 x(3) z(3);]');
                    bb(1)=(z(2)-z(3))/De;  cc(1)=(x(3)-x(2))/De;  aa(1)=(x(2)*z(3)-x(3)*z(2))/De;
                    bb(2)=(z(3)-z(1))/De;  cc(2)=(x(1)-x(3))/De;  aa(2)=(x(3)*z(1)-x(1)*z(3))/De;
                    bb(3)=(z(1)-z(2))/De;  cc(3)=(x(2)-x(1))/De;  aa(3)=(x(1)*z(2)-x(2)*z(1))/De;
                    %----------------------------------------------------------
                    zeta(1)=aa(1)+bb(1)*xx(ii,jj)+cc(1)*zz(ii,jj);
                    zeta(2)=aa(2)+bb(2)*xx(ii,jj)+cc(2)*zz(ii,jj);
                    zeta(3)=aa(3)+bb(3)*xx(ii,jj)+cc(3)*zz(ii,jj);
                end
            end,edges=[TModel.Edges(element.Edges)];ll=[edges.Length];
            wwx(1)=ll(1)*(zeta(1)*bb(2)-zeta(2)*bb(1));wwz(1)=ll(1)*(zeta(1)*cc(2)-zeta(2)*cc(1));wwx(2)=ll(2)*(zeta(2)*bb(3)-zeta(3)*bb(2));
            wwz(2)=ll(2)*(zeta(2)*cc(3)-zeta(3)*cc(2));wwx(3)=ll(3)*(zeta(3)*bb(1)-zeta(1)*bb(3));wwz(3)=ll(3)*(zeta(3)*cc(1)-zeta(1)*cc(3));
            for kk=1:3
                Ex(ii,jj)=Ex(ii,jj)+element.EdgeSigns(kk)*wwx(kk)*EdgeVector_E(edges(kk).Index2D);Ez(ii,jj)=Ez(ii,jj)+element.EdgeSigns(kk)*wwz(kk)*EdgeVector_E(edges(kk).Index2D);
                Ey(ii,jj)=Ey(ii,jj)+zeta(kk)*VertexVector_E(vertices(kk).Index2D);
            end
        end
    end
    [cmap]=buildcmap('cbkry');figure;title("Real Fields");
    subplot(1,3,1);pcolor(yy,zz,real(Ex));shading interp;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(yy,zz,real(Ey));shading interp;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(yy,zz,real(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
    [cmap]=buildcmap('cbkry');figure;title("Imaginary Fields");
    subplot(1,3,1);pcolor(yy,zz,imag(Ex));shading interp;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(yy,zz,imag(Ey));shading interp;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(yy,zz,imag(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
    [cmap]=buildcmap('cbkry');figure;title("Absolute Fields");
    subplot(1,3,1);pcolor(yy,zz,abs(Ex));shading interp;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(yy,zz,abs(Ey));shading interp;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(yy,zz,abs(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;

end
function [] = Plot2D_E_TFLF_z(TModel,BoundaryExcitation,Ex,Ey,Ez,VertexVector_E,EdgeVector_E,xx,yy),ci=0;
    elements=[TModel.Facets(BoundaryExcitation.Facets)];BaryCenters=[elements.Barycenter];xbc=BaryCenters(1,:);ybc=BaryCenters(2,:);
    for ii=1:size(xx,1)
        for jj=1:size(yy,2),ic=0;ie=0;dist=(xx(ii,jj)-xbc).^2 +(yy(ii,jj)-ybc).^2;[~,Is]=sort(dist);
            while ie==0,ic=ic+1;element=elements(Is(ic));vertices=[TModel.Vertices(element.Vertices)];x=[vertices.X];y=[vertices.Y];
                De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');
                bb(1)=(y(2)-y(3))/De;  cc(1)=(x(3)-x(2))/De;  aa(1)=(x(2)*y(3)-x(3)*y(2))/De;
                bb(2)=(y(3)-y(1))/De;  cc(2)=(x(1)-x(3))/De;  aa(2)=(x(3)*y(1)-x(1)*y(3))/De;
                bb(3)=(y(1)-y(2))/De;  cc(3)=(x(2)-x(1))/De;  aa(3)=(x(1)*y(2)-x(2)*y(1))/De;
                %----------------------------------------------------------
                zeta(1)=aa(1)+bb(1)*xx(ii,jj)+cc(1)*yy(ii,jj);
                zeta(2)=aa(2)+bb(2)*xx(ii,jj)+cc(2)*yy(ii,jj);
                zeta(3)=aa(3)+bb(3)*xx(ii,jj)+cc(3)*yy(ii,jj);
                if (zeta(1)>0 || abs(zeta(1))<1000*eps) && (zeta(2)>0 || abs(zeta(2))<1000*eps) && (zeta(3)>0 || abs(zeta(3))<1000*eps),ie=Is(ic);end
                if(ic==numel(Is)),ie=Is(1);ci=ci+1;%disp(ci);
                    element=elements(Is(1));vertices=[TModel.Vertices(element.Vertices)];x=[vertices.X];y=[vertices.Y];
                    De=det([1 x(1) y(1);1 x(2) y(2);1 x(3) y(3);]');
                    bb(1)=(y(2)-y(3))/De;  cc(1)=(x(3)-x(2))/De;  aa(1)=(x(2)*y(3)-x(3)*y(2))/De;
                    bb(2)=(y(3)-y(1))/De;  cc(2)=(x(1)-x(3))/De;  aa(2)=(x(3)*y(1)-x(1)*y(3))/De;
                    bb(3)=(y(1)-y(2))/De;  cc(3)=(x(2)-x(1))/De;  aa(3)=(x(1)*y(2)-x(2)*y(1))/De;
                    %----------------------------------------------------------
                    zeta(1)=aa(1)+bb(1)*xx(ii,jj)+cc(1)*yy(ii,jj);
                    zeta(2)=aa(2)+bb(2)*xx(ii,jj)+cc(2)*yy(ii,jj);
                    zeta(3)=aa(3)+bb(3)*xx(ii,jj)+cc(3)*yy(ii,jj);
                end
            end,edges=[TModel.Edges(element.Edges)];ll=[edges.Length];
            wwx(1)=ll(1)*(zeta(1)*bb(2)-zeta(2)*bb(1));wwy(1)=ll(1)*(zeta(1)*cc(2)-zeta(2)*cc(1));wwx(2)=ll(2)*(zeta(2)*bb(3)-zeta(3)*bb(2));
            wwy(2)=ll(2)*(zeta(2)*cc(3)-zeta(3)*cc(2));wwx(3)=ll(3)*(zeta(3)*bb(1)-zeta(1)*bb(3));wwy(3)=ll(3)*(zeta(3)*cc(1)-zeta(1)*cc(3));
            for kk=1:3
                Ex(ii,jj)=Ex(ii,jj)+element.EdgeSigns(kk)*wwx(kk)*EdgeVector_E(edges(kk).Index2D);Ey(ii,jj)=Ey(ii,jj)+element.EdgeSigns(kk)*wwy(kk)*EdgeVector_E(edges(kk).Index2D);
                Ez(ii,jj)=Ez(ii,jj)+zeta(kk)*VertexVector_E(vertices(kk).Index2D);
            end
        end
    end
    [cmap]=buildcmap('cbkry');figure;title("Real Fields");
    subplot(1,3,1);pcolor(yy,zz,real(Ex));shading interp;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(yy,zz,real(Ey));shading interp;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(yy,zz,real(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
    [cmap]=buildcmap('cbkry');figure;title("Imaginary Fields");
    subplot(1,3,1);pcolor(yy,zz,imag(Ex));shading interp;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(yy,zz,imag(Ey));shading interp;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(yy,zz,imag(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
    [cmap]=buildcmap('cbkry');figure;title("Absolute Fields");
    subplot(1,3,1);pcolor(yy,zz,abs(Ex));shading interp;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,2);pcolor(yy,zz,abs(Ey));shading interp;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
    subplot(1,3,3);pcolor(yy,zz,abs(Ez));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
end