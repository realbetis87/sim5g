%--------------------------------------------------------------------------
%{
                    Plot E-B Excitation Formulation Fields

       1. X          : Uknown Field to be plotted (1 x Total Number Of  Uknowns)
       2. Position   : Plane Coordinates - for (-Ax) plane
                               Position=[-Ax;Nan;Nan;]
       3. Points     : [NumberOfPoints 1 Axis, NumberOfPOints 2 Axis]

%}
%--------------------------------------------------------------------------
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz] = ReturnFields_Excitation(TModel,X,Position,Points),ElectromagneticConstants;
    NumberOfFacets=numel(TModel.Facets);NumberOfEdges=numel(TModel.Edges);FacetVector=zeros(NumberOfFacets,1);EdgeVector=zeros(NumberOfEdges,1);
    for ii=1:NumberOfFacets,facet=TModel.Facets(ii);if(facet.UknownIndex~=0),FacetVector(ii)=-1i*X(abs(facet.UknownIndex))/c0;end,end
    for ii=1:NumberOfEdges,edge=TModel.Edges(ii);if(edge.UknownIndex~=0),EdgeVector(ii)=X(abs(edge.UknownIndex));end,end
    Plane=find(~isnan(Position));N1=Points(1);N2=Points(2);empty=zeros(N1,N2);Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
    vertices=[TModel.Vertices];xs=[vertices.X];ys=[vertices.Y];zs=[vertices.Z];
    Xmin=min(xs);Ymin=min(ys);Zmin=min(zs);Xmax=max(xs);Ymax=max(ys);Zmax=max(zs);
    switch Plane
        case 1,xx=Position(1);y=linspace(Ymin-Ymin*1e-3,Ymax-Ymax*1e-3,N1);z=linspace(Zmin-Zmin*1e-3,Zmax-Zmax*1e-3,N2);[yy,zz]=meshgrid(y,z);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotXPlane_Field(TModel,EdgeVector,FacetVector,xx,yy,zz);
        case 2,x=linspace(Xmin-Xmin*1e-3,Xmax-Xmax*1e-3,N1);yy=Position(2);z=linspace(Zmin-Zmin*1e-3,Zmax-Zmax*1e-3,N2);[xx,zz]=meshgrid(x,z);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotYPlane_Field(TModel,EdgeVector,FacetVector,xx,yy,zz);
        case 3,x=linspace(Xmin-Xmin*1e-3,Xmax-Xmax*1e-3,N1);y=linspace(Ymin-Ymin*1e-3,Ymax-Ymax*1e-3,N2);zz=Position(3);[xx,yy]=meshgrid(x,y);[Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotZPlane_Field(TModel,EdgeVector,FacetVector,xx,yy,zz);
     end
end
%--------------------------------------------------------------------------
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotXPlane_Field(TModel,EdgeVector,FacetVector,xx,yy,zz),error=1000*eps;empty=zeros(size(yy,1),size(yy,2));Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
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

            dfi=zeros(4,1);dei=zeros(6,1);
            for mm=1:6,if(edges(mm).UknownIndex==0),dei(mm)=1;elseif(sign(edges(mm).UknownIndex)<0),dei(mm)=-1;else,dei(mm)=1;end,end
            for mm=1:4,if(facets(mm).UknownIndex==0),dfi(mm)=1;elseif(sign(facets(mm).UknownIndex)<0),dfi(mm)=-1;else,dfi(mm)=1;end,end
            for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index);
                                                  By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index);
                                                  Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index);
            end
            for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index);
                                                 Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index);
                                                 Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index);
            end
         end
    end
end
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotYPlane_Field(TModel,EdgeVector,FacetVector,xx,yy,zz),error=1000*eps;empty=zeros(size(xx,1),size(xx,2));Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
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

            dfi=zeros(4,1);dei=zeros(6,1);
            for mm=1:6,if(edges(mm).UknownIndex==0),dei(mm)=1;elseif(sign(edges(mm).UknownIndex)<0),dei(mm)=-1;else,dei(mm)=1;end,end
            for mm=1:4,if(facets(mm).UknownIndex==0),dfi(mm)=1;elseif(sign(facets(mm).UknownIndex)<0),dfi(mm)=-1;else,dfi(mm)=1;end,end
            for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index);
                                                  By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index);
                                                  Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index);
            end
            for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index);
                                                 Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index);
                                                 Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index);
            end
         end
    end
end
function [Ex,Ey,Ez,Bx,By,Bz,xx,yy,zz]=PlotZPlane_Field(TModel,EdgeVector,FacetVector,xx,yy,zz),error=1000*eps;empty=zeros(size(xx,1),size(xx,2));Ex=empty;Ey=empty;Ez=empty;Bx=empty;By=empty;Bz=empty;
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

            dfi=zeros(4,1);dei=zeros(6,1);
            for mm=1:6,if(edges(mm).UknownIndex==0),dei(mm)=1;elseif(sign(edges(mm).UknownIndex)<0),dei(mm)=-1;else,dei(mm)=1;end,end
            for mm=1:4,if(facets(mm).UknownIndex==0),dfi(mm)=1;elseif(sign(facets(mm).UknownIndex)<0),dfi(mm)=-1;else,dfi(mm)=1;end,end
            for kk=1:4,sgn=element.FacetSigns(kk);Bx(ii,jj)=Bx(ii,jj)+wfx(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index);
                                                  By(ii,jj)=By(ii,jj)+wfy(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index);
                                                  Bz(ii,jj)=Bz(ii,jj)+wfz(kk)*dfi(kk)*sgn*FacetVector(facets(kk).Index);
            end
            for kk=1:6,sgn=element.EdgeSigns(kk);Ex(ii,jj)=Ex(ii,jj)+wx(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index);
                                                 Ey(ii,jj)=Ey(ii,jj)+wy(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index);
                                                 Ez(ii,jj)=Ez(ii,jj)+wz(kk)*dei(kk)*sgn*EdgeVector(edges(kk).Index);
            end
         end
    end
end
