%Olusturulmus olan son nokta bulutunun girdileri
A_son=300;
B_son=round(A_son/3);
C_son=B_son;
D_son=A_son-B_son-C_son;

son=zeros(3,A_son);
son(1:2,1:B_son)=rand(2,B_son);
son([1,3],(B_son+1):(B_son+C_son))=rand(2,C_son);
son(2,(B_son+1):(B_son+C_son))=1;
son(2:3,(B_son+C_son+1):A_son)=rand(2,D_son);
son(1,(B_son+C_son+1):A_son)=1;
%Alinan verinin nokta bulut girdileri
A_veri=100;
B_veri=round(A_veri/3);
C_veri=B_veri;
D_veri=A_veri-B_veri-C_veri;

veri=zeros(3,A_veri);
veri(1:2,1:B_veri)=rand(2,B_veri);
veri([1,3],(B_veri+1):(B_veri+C_veri))=rand(2,C_veri);
veri(2,(B_veri+1):(B_veri+C_veri))=1;
veri(2:3,(B_veri+C_veri+1):A_veri)=rand(2,D_veri);
veri(1,(B_veri+C_veri+1):A_veri)=1;

%veri noktalarinin baslangic halleri

x_konum=2*(2*rand-1); 
y_konum=5*(2*rand-1); 
z_konum=2*(2*rand-1);
X=[1 0 0;0 cos(x_konum) -sin(x_konum);0 sin(x_konum) cos(x_konum)];
Y=[cos(y_konum) 0 sin(y_konum);0 1 0;-sin(y_konum) 0 cos(y_konum)];
Z=[cos(z_konum) -sin(z_konum) 0;sin(z_konum) cos(z_konum) 0;0 0 1];

R=Z*Y*X;
veri=R*veri;
veri(1,:)=veri(1,:)+0.2*randn;
veri(2,:)=veri(2,:)+0.2*randn;
veri(3,:)=veri(3,:)+0.2*randn;

figure(1)
plot3(son(1,:),son(2,:),son(3,:),'b.',veri(1,:),veri(2,:),veri(3,:),'r.'), hold on, axis equal
plot3([1 1 0],[0 1 1],[0 0 0],'r-',[1 1],[1 1],[0 1],'r-','LineWidth',2)
title('Alinan veriler (kirmizi) ve modellenmis noktalar (mavi)')

[RotMat,TransVec,dataOut]=icp(son,veri);

figure(2)
plot3(son(1,:),son(2,:),son(3,:),'b.',dataOut(1,:),dataOut(2,:),dataOut(3,:),'y.'), hold on, axis equal
plot3([1 1 0],[0 1 1],[0 0 0],'r-',[1 1],[1 1],[0 1],'r-','LineWidth',2)
title('Donusturulmus veriler (sari) ve modellenmis noktalar (mavi)')

