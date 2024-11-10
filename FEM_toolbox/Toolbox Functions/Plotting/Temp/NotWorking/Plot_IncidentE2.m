function []= Plot_IncidentE2(Vertices,Elements,Field,xx,yy,zz,Element_Edges,Edge_id,known_E)
error=1000*eps;NumberOfElements=size(Elements,2);
NumberOfEdges=size(Edge_id,1);
St=Field;

for in=1:NumberOfEdges
%    St(in)=Field(known_E(in));
    if(Edge_id(in)~=2)
        St(in)=Field(known_E(in));
    end
end
%}
Ny=size(yy,2);Nz=size(zz,2);
Ex=zeros(Ny,Nz);Ey=zeros(Ny,Nz);Ez=zeros(Ny,Nz);
xg=(Vertices(1,Elements(1,1:NumberOfElements))+Vertices(1,Elements(2,1:NumberOfElements))+Vertices(1,Elements(3,1:NumberOfElements))+Vertices(1,Elements(4,1:NumberOfElements)))/4;
yg=(Vertices(2,Elements(1,1:NumberOfElements))+Vertices(2,Elements(2,1:NumberOfElements))+Vertices(2,Elements(3,1:NumberOfElements))+Vertices(2,Elements(4,1:NumberOfElements)))/4;
zg=(Vertices(3,Elements(1,1:NumberOfElements))+Vertices(3,Elements(2,1:NumberOfElements))+Vertices(3,Elements(3,1:NumberOfElements))+Vertices(3,Elements(4,1:NumberOfElements)))/4;
for ii=1:Ny
    for jj=1:Nz
        aux=0;ie=0;dist=(xx-xg).^2 +(yy(ii)-yg).^2+(zz(jj)-zg).^2;[distsort,Is]=sort(dist);
        while ie==0
             aux=aux+1;nodes(1:4)=Elements(1:4,Is(aux));
             x(1:4)=Vertices(1,nodes(1:4));y(1:4)=Vertices(2,nodes(1:4));zd(1:4)=Vertices(3,nodes(1:4));
             De=det([1 x(1) y(1) zd(1);1 x(2) y(2) zd(2);1 x(3) y(3) zd(3);1 x(4) y(4) zd(4);]);
             a(1)=det ([1 x(1) y(1) zd(1); 0 x(2) y(2) zd(2); 0 x(3) y(3) zd(3); 0 x(4) y(4) zd(4)]) / De;
             a(2)=det ([0 x(1) y(1) zd(1); 1 x(2) y(2) zd(2); 0 x(3) y(3) zd(3); 0 x(4) y(4) zd(4)]) / De;
             a(3)=det ([0 x(1) y(1) zd(1); 0 x(2) y(2) zd(2); 1 x(3) y(3) zd(3); 0 x(4) y(4) zd(4)]) / De;
             a(4)=det ([0 x(1) y(1) zd(1); 0 x(2) y(2) zd(2); 0 x(3) y(3) zd(3); 1 x(4) y(4) zd(4)]) / De;
             b(1)=det ([1 1 y(1) zd(1); 1 0 y(2) zd(2); 1 0 y(3) zd(3); 1 0 y(4) zd(4)]) / De;         
             b(2)=det ([1 0 y(1) zd(1); 1 1 y(2) zd(2); 1 0 y(3) zd(3); 1 0 y(4) zd(4)]) / De;          
             b(3)=det ([1 0 y(1) zd(1); 1 0 y(2) zd(2); 1 1 y(3) zd(3); 1 0 y(4) zd(4)]) / De;   
             b(4)=det ([1 0 y(1) zd(1); 1 0 y(2) zd(2); 1 0 y(3) zd(3); 1 1 y(4) zd(4)]) / De;  
             c(1)=det ([1 x(1) 1 zd(1); 1 x(2) 0 zd(2); 1 x(3) 0 zd(3); 1 x(4) 0 zd(4)]) / De;         
             c(2)=det ([1 x(1) 0 zd(1); 1 x(2) 1 zd(2); 1 x(3) 0 zd(3); 1 x(4) 0 zd(4)]) / De;         
             c(3)=det ([1 x(1) 0 zd(1); 1 x(2) 0 zd(2); 1 x(3) 1 zd(3); 1 x(4) 0 zd(4)]) / De;
             c(4)=det ([1 x(1) 0 zd(1); 1 x(2) 0 zd(2); 1 x(3) 0 zd(3); 1 x(4) 1 zd(4)]) / De;
             d(1)=det ([1 x(1) y(1) 1; 1 x(2) y(2) 0; 1 x(3) y(3) 0; 1 x(4) y(4) 0]) / De;
             d(2)=det ([1 x(1) y(1) 0; 1 x(2) y(2) 1; 1 x(3) y(3) 0; 1 x(4) y(4) 0]) / De;
             d(3)=det ([1 x(1) y(1) 0; 1 x(2) y(2) 0; 1 x(3) y(3) 1; 1 x(4) y(4) 0]) / De;
             d(4)=det ([1 x(1) y(1) 0; 1 x(2) y(2) 0; 1 x(3) y(3) 0; 1 x(4) y(4) 1]) / De;
             z(1)=a(1)+b(1)*xx+c(1)*yy(ii)+d(1)*zz(jj);z(2)=a(2)+b(2)*xx+c(2)*yy(ii)+d(2)*zz(jj);
             z(3)=a(3)+b(3)*xx+c(3)*yy(ii)+d(3)*zz(jj);z(4)=a(4)+b(4)*xx+c(4)*yy(ii)+d(4)*zz(jj);
             ll(1)=sqrt((x(1)-x(2))^2 +(y(1)-y(2))^2 +(zd(1)-zd(2))^2);
             ll(2)=sqrt((x(1)-x(3))^2 +(y(1)-y(3))^2 +(zd(1)-zd(3))^2);
             ll(3)=sqrt((x(1)-x(4))^2 +(y(1)-y(4))^2 +(zd(1)-zd(4))^2);
             ll(4)=sqrt((x(2)-x(3))^2 +(y(2)-y(3))^2 +(zd(2)-zd(3))^2);
             ll(5)=sqrt((x(2)-x(4))^2 +(y(2)-y(4))^2 +(zd(2)-zd(4))^2);
             ll(6)=sqrt((x(3)-x(4))^2 +(y(3)-y(4))^2 +(zd(3)-zd(4))^2);
             ll=ones(6,1);
             if (z(1)>0 || abs(z(1))<=error) && (z(2)>0 || abs(z(2))<=error) && (z(3)>0 || abs(z(3))<=error)&& (z(4)>0 || abs(z(4))<=error)
              ie=Is(aux);
             end
        end
              wx(1)=z(1)*b(2)-z(2)*b(1); wx(2)=z(1)*b(3)-z(3)*b(1);wx(3)=z(1)*b(4)-z(4)*b(1);wx(4)=z(2)*b(3)-z(3)*b(2);wx(5)=z(2)*b(4)-z(4)*b(2);wx(6)=z(3)*b(4)-z(4)*b(3);
              wy(1)=z(1)*c(2)-z(2)*c(1); wy(2)=z(1)*c(3)-z(3)*c(1);wy(3)=z(1)*c(4)-z(4)*c(1);wy(4)=z(2)*c(3)-z(3)*c(2);wy(5)=z(2)*c(4)-z(4)*c(2);wy(6)=z(3)*c(4)-z(4)*c(3);
              wz(1)=z(1)*d(2)-z(2)*d(1); wz(2)=z(1)*d(3)-z(3)*d(1);wz(3)=z(1)*d(4)-z(4)*d(1);wz(4)=z(2)*d(3)-z(3)*d(2);wz(5)=z(2)*d(4)-z(4)*d(2);wz(6)=z(3)*d(4)-z(4)*d(3);
              wx(1)=wx(1)*ll(1);wy(1)=wy(1)*ll(1);wz(1)=wz(1)*ll(1);    
              wx(2)=wx(2)*ll(2);wy(2)=wy(2)*ll(2);wz(2)=wz(2)*ll(2);
              wx(3)=wx(3)*ll(3);wy(3)=wy(3)*ll(3);wz(3)=wz(3)*ll(3); 
              wx(4)=wx(4)*ll(4);wy(4)=wy(4)*ll(4);wz(4)=wz(4)*ll(4); 
              wx(5)=wx(5)*ll(5);wy(5)=wy(5)*ll(5);wz(5)=wz(5)*ll(5); 
              wx(6)=wx(6)*ll(6);wy(6)=wy(6)*ll(6);wz(6)=wz(6)*ll(6); 
              edge(1:6)=Element_Edges(1:6,ie);
              for kk=1:6
                   %if (Edge_id(abs(edge(kk)))==1)
                    sgn=sign(edge(kk));
                    Ex(ii,jj)=Ex(ii,jj)+wx(kk)*sgn*St(abs(edge(kk)));
                    Ey(ii,jj)=Ey(ii,jj)+wy(kk)*sgn*St(abs(edge(kk)));
                    Ez(ii,jj)=Ez(ii,jj)+wz(kk)*sgn*St(abs(edge(kk)));
                   %end
              end
    end
end
[cmap]=buildcmap('cbkry');figure;
subplot(1,3,1);pcolor(real(Ex'));shading interp;hold on;colormap(cmap);title('x Component Electric Field');colorbar;axis tight;axis equal;
subplot(1,3,2);pcolor(real(Ey'));shading interp;hold on;colormap(cmap);title('y Component Electric Field');colorbar;axis tight;axis equal;
subplot(1,3,3);pcolor(real(Ez'));shading interp;colormap(cmap);title('z Component Electric Field');colorbar;axis tight;axis equal;
