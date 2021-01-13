
m = 80; %genislik
n = m^2; % nokta sayisi
[X,Y] = meshgrid(linspace(-2,2,m), linspace(-2,2,m));
X = reshape(X,1,[]);
Y = reshape(Y,1,[]);
Z = sin(X).*cos(Y);
% Veri nokta matrisi
Veri = [X; Y; Z];

% cevrim
Cevrim_x = 0.5;
Cevrim_y = -0.3;
Cevrim_z = 0.2;

% cevrim vektor
T = [Cevrim_x; Cevrim_y; Cevrim_z];

% rotasyon miktari
Rotasyon_x = 0.3;
Rotasyon_y = -0.2;
Rotasyon_z = 0.05;

Rotas_x = [1 0 0;
      0 cos(Rotasyon_x) -sin(Rotasyon_x);
      0 sin(Rotasyon_x) cos(Rotasyon_x)];
  
Rotas_y = [cos(Rotasyon_y) 0 sin(Rotasyon_y);
      0 1 0;
      -sin(Rotasyon_y) 0 cos(Rotasyon_y)];
  
Rotas_z = [cos(Rotasyon_z) -sin(Rotasyon_z) 0;
      sin(Rotasyon_z) cos(Rotasyon_z) 0;
      0 0 1];

% Rotasyon matrisi
R = Rotas_x*Rotas_y*Rotas_z;

% veri matris donusturme ve gurultu ekleme
Model = R * Veri + repmat(T, 1, n);

% model ve veriye gurultu
rng(100);
Model = Model + 0.01*randn(3,n);
Veri = Veri + 0.01*randn(3,n);

[Ricp Ticp ER t] = icp(Model, Veri, 15);

Dicp = Ricp * Veri + repmat(Ticp, 1, n);

figure;
subplot(2,2,1);
plot3(Model(1,:),Model(2,:),Model(3,:),'bo',Veri(1,:),Veri(2,:),Veri(3,:),'y.');
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
title('sari: z=sin(x)*cos(y), mavi: donusturulmus nokta bulutu');

subplot(2,2,2);
plot3(Model(1,:),Model(2,:),Model(3,:),'bo',Dicp(1,:),Dicp(2,:),Dicp(3,:),'y.');
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
title('ICP');

subplot(2,2,[3 4]);
plot(0:15,ER,'--x');
xlabel('iterasyon sayisi');
ylabel('hata');
legend('ICP algoritmasi');
title(['Gecen vakit: ' num2str(t(end),2) ' s']);

[Ricp Ticp ER t] = icp(Model, Veri, 15, 'Matching', 'kDtree', 'Extrapolation', true);

Dicp = Ricp * Veri + repmat(Ticp, 1, n);


figure;
subplot(2,2,1);
plot3(Model(1,:),Model(2,:),Model(3,:),'bo',Veri(1,:),Veri(2,:),Veri(3,:),'y.');
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
title('sari: z=sin(x)*cos(y), mavi: donusturulmus nokta bulutu');


subplot(2,2,2);
plot3(Model(1,:),Model(2,:),Model(3,:),'bo',Dicp(1,:),Dicp(2,:),Dicp(3,:),'y.');
axis equal;
xlabel('x'); ylabel('y'); zlabel('z');
title('ICP');

subplot(2,2,[3 4]);
plot(0:15,ER,'--x');
xlabel('iterasyon sayisi');
ylabel('hata');
legend('k-d agaci');
title(['Gecen sure: ' num2str(t(end),2) ' s']);
