function [] = Plot2DElectricField(Boundary,Vertices2D,Edges2D,Elements2D,Input),VertexField=zeros(numel(Vertices2D),1);EdgeField=zeros(numel(Edges2D),1);
%==========================================================================
%{
        Plot Function for 2D Electric Field Input 

%}
%==========================================================================
    for ii=1:numel(Vertices2D),vertex=Vertices2D(ii);if(vertex.Id~='e'),VertexField(ii)=Input(vertex.EIndex);end,end
    for ii=1:numel(Edges2D),edge=Edges2D(ii);if(edge.Id~='e'),EdgeField(ii)=Input(edge.EIndex);end,end
    if(Boundary.Geometry.Axis == 'x'),PlotEFieldAtX(Boundary,Vertices2D,Edges2D,Elements2D,VertexField,EdgeField);
    elseif(Boundary.Geometry.Axis=='y'),PlotEFieldAtY(Boundary,Vertices2D,Edges2D,Elements2D,VertexField,EdgeField);
    elseif(Boundary.Geometry.Axis=='z'),PlotEFieldAtZ(Boundary,Vertices2D,Edges2D,Elements2D,VertexField,EdgeField);
    end
end

function [] = PlotEFieldAtX(Boundary,Vertices2D,Edges2D,Element2D,VertexField,EdgeField),Ny=100;Nz=100;geometry=Boundary.Geometry;[yy,zz]=ReturnPoints(geometry,Ny,Nz);Ex=zeros(Ny,Nz);Ey=zeros(Ny,Nz);Ez=zeros(Ny,Nz);
    %Find Centroids Cx,Cy,Cz
    for ii=1:Ny
        for jj=1:Nz
            
        end
    end
   subplot(2,3,1);pcolor(yy,zz,real(Ex));shading interp;hold on;colormap(cmap);title('x Component Real Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,2);pcolor(yy,zz,real(Ey));shading interp;hold on;colormap(cmap);title('y Component Real Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,3);pcolor(yy,zz,real(Ez));shading interp;hold on;colormap(cmap);title('z Component Real Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,4);pcolor(yy,zz,imag(Ex));shading interp;hold on;colormap(cmap);title('x Component Imaginary Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,5);pcolor(yy,zz,imag(Ey));shading interp;hold on;colormap(cmap);title('y Component Imaginary Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,6);pcolor(yy,zz,imag(Ez));shading interp;hold on;colormap(cmap);title('z Component Imaginary Electric Field');colorbar;axis tight;axis equal;
end
function [] = PlotEFieldAtY(Boundary,Vertices2D,Edges2D,Element2D,VertexField,EdgeField),Nx=100;Nz=100;geometry=Boundary.Geometry;[xx,zz]=ReturnPoints(geometry,Nx,Nz);Ex=zeros(Nx,Nz);Ey=zeros(Nx,Nz);Ez=zeros(Nx,Nz);
    for ii=1:Nx
        for jj=1:Nz
        end
    end
   subplot(2,3,1);pcolor(xx,zz,real(Ex));shading interp;hold on;colormap(cmap);title('x Component Real Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,2);pcolor(xx,zz,real(Ey));shading interp;hold on;colormap(cmap);title('y Component Real Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,3);pcolor(xx,zz,real(Ez));shading interp;hold on;colormap(cmap);title('z Component Real Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,4);pcolor(xx,zz,imag(Ex));shading interp;hold on;colormap(cmap);title('x Component Imaginary Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,5);pcolor(xx,zz,imag(Ey));shading interp;hold on;colormap(cmap);title('y Component Imaginary Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,6);pcolor(xx,zz,imag(Ez));shading interp;hold on;colormap(cmap);title('z Component Imaginary Electric Field');colorbar;axis tight;axis equal;
end
function [] = PlotEFieldAtZ(Boundary,Vertices2D,Edges2D,Element2D,VertexField,EdgeField),Nx=100;Ny=100;geometry=Boundary.Geometry;[xx,yy]=ReturnPoints(geometry,Nx,Ny);Ex=zeros(Nx,Ny);Ey=zeros(Nx,Ny);Ez=zeros(Nx,Ny);
    for ii=1:Nx
        for jj=1:Ny
        end
    end
   subplot(2,3,1);pcolor(xx,yy,real(Ex));shading interp;hold on;colormap(cmap);title('x Component Real Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,2);pcolor(xx,yy,real(Ey));shading interp;hold on;colormap(cmap);title('y Component Real Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,3);pcolor(xx,yy,real(Ez));shading interp;hold on;colormap(cmap);title('z Component Real Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,4);pcolor(xx,yy,imag(Ex));shading interp;hold on;colormap(cmap);title('x Component Imaginary Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,5);pcolor(xx,yy,imag(Ey));shading interp;hold on;colormap(cmap);title('y Component Imaginary Electric Field');colorbar;axis tight;axis equal;
   subplot(2,3,6);pcolor(xx,yy,imag(Ez));shading interp;hold on;colormap(cmap);title('z Component Imaginary Electric Field');colorbar;axis tight;axis equal;
end


function [HorizontalPoints,VerticalPoints] = ReturnPoints(Geometry,NumberOfPointsHorizontal,NumberOfPointsVertical)
    if(isa(Geometry,"Plane"))
        switch Geometry.Axis
            case 'x',Ymin=min(min(Geometry.YCorners));Ymax=max(max(Geometry.YCorners));Zmin=min(min(Geometry.ZCorners));Zmax=max(max(Geometry.ZCorners));HorizontalDimension=Ymax-Ymin;VerticalDimension=Zmax-Zmin;
                     HorizontalPoints=linspace(Ymin,HorizontalDimension/NumberOfPointsHorizontal,Ymax);VerticalPoints=linspace(Zmin,VerticalDimension/NumberOfPointsVertical,Zmax);
            case 'y',Xmin=min(min(Geometry.XCorners));Xmax=max(max(Geometry.XCorners));Zmin=min(min(Geometry.ZCorners));Zmax=max(max(Geometry.ZCorners));HorizontalDimension=Xmax-Xmin;VerticalDimension=Zmax-Zmin;
                     HorizontalPoints=linspace(Xmin,HorizontalDimension/NumberOfPointsHorizontal,Xmax);VerticalPoints=linspace(Zmin,VerticalDimension/NumberOfPointsVertical,Zmax);
            case 'z',Ymin=min(min(Geometry.YCorners));Ymax=max(max(Geometry.YCorners));Xmin=min(min(Geometry.XCorners));Xmax=max(max(Geometry.XCorners));HorizontalDimension=Xmax-Xmin;VerticalDimension=Ymax-Ymin;
                     HorizontalPoints=linspace(Xmin,HorizontalDimension/NumberOfPointsHorizontal,Xmax);VerticalPoints=linspace(Ymin,VerticalDimension/NumberOfPointsVertical,Ymax);
        end
    elseif(isa(Geometry,"Circle2D")),phi=linspace(0,2*pi,NumberOfPointsVer);rad=linspace((1/1000)*geometry.Radius,(geometry.Radius)-(1/200)*geometry.Radius,NumberOfPointsHor);[R,Phi]=meshgrid(rad,phi);xx=R.*cos(Phi);yy=R.*sin(Phi);
    end
end
