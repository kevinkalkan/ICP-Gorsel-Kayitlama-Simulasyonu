function [TR, TT, Hata, t_vektor] = icp(q,p,varargin)
% [TR, TT] = icp(q,p)   TR rotasyon matrisini döndürür 
% TT çevrim vektörünü verir (TR * p + TT ) ile q arasindaki mesafeyi kisaltir 
% p 3xM matris, q 3xN matris
% [TR, TT] = icp(q,p,k)   algoritman?n k iterasyonunda yapilmasi
% [TR, TT, ER] = icp(q,p,k)  k iterasyonda RMS hatasi gönderir
% [TR, TT, ER, t] = icp(q,p,k)   k iterasyondaki t süresini belirtir

inp = inputParser;
inp.addRequired('q', @(x)isreal(x) && size(x,1) == 3);
inp.addRequired('p', @(x)isreal(x) && size(x,1) == 3);

inp.addOptional('iter', 10, @(x)x > 0 && x < 10^5);

inp.addParamValue('Boundary', [], @(x)size(x,1) == 1);
inp.addParamValue('EdgeRejection', false, @(x)islogical(x)); 
inp.addParamValue('Extrapolation', false, @(x)islogical(x));
validMatching = {'bruteForce','Delaunay','kDtree'};
inp.addParamValue('Matching', 'bruteForce', @(x)any(strcmpi(x,validMatching)));
validMinimize = {'point','plane','lmapoint'};
inp.addParamValue('Minimize', 'point', @(x)any(strcmpi(x,validMinimize)));
inp.addParamValue('Normals', [], @(x)isreal(x) && size(x,1) == 3);
inp.addParamValue('NormalsData', [], @(x)isreal(x) && size(x,1) == 3);
inp.addParamValue('ReturnAll', false, @(x)islogical(x));
inp.addParamValue('Triangulation', [], @(x)isreal(x) && size(x,2) == 3);
inp.addParamValue('Verbose', false, @(x)islogical(x));
inp.addParamValue('Weight', @(x)ones(1,length(x)), @(x)isa(x,'function_handle'));
inp.addParamValue('WorstRejection', 0, @(x)isscalar(x) && x > 0 && x < 1);

inp.parse(q,p,varargin{:});
arg = inp.Results;
clear('inp');

% RMS hatasi icin vektore yer ayrilir
t_vektor = zeros(arg.iter+1,1); 

% Start timer
tic;

Np = size(p,2);

% donusturulmus veri nokta bulutu
nokta = p;

Hata = zeros(arg.iter+1,1); 

% matris ve vektor
T = zeros(3,1);
R = eye(3,3);

% rotasyon matrisi ve cevrim vektoru
TT = zeros(3,1, arg.iter+1);
TR = repmat(eye(3,3), [1,1, arg.iter+1]);
    
if (strcmp(arg.Minimize, 'plane') && isempty(arg.Normals))
    arg.Normals = lsqnormest(q,4); 
end

if strcmp(arg.Matching, 'Delaunay')
    DT = DelaunayTri(transpose(q));
end

if strcmp(arg.Matching, 'kDtree')
    kdOBJ = KDTreeSearcher(transpose(q));
end

if arg.EdgeRejection
    if isempty(arg.Boundary)
        bdr = find_bound(q, arg.Triangulation);
    else
        bdr = arg.Boundary;
    end
end

if arg.Extrapolation
    qq = [ones(1,arg.iter+1);zeros(6,arg.iter+1)];      
    dq = zeros(7,arg.iter+1);
    theta = zeros(1,arg.iter+1);
end

t_vektor(1) = toc;

% iterasyon
for k=1:arg.iter
    switch arg.Matching
        case 'bruteForce'
            [match dagilim_min] = match_bruteForce(q,nokta);
        case 'Delaunay'
            [match dagilim_min] = match_Delaunay(q,nokta,DT);
        case 'kDtree'
            [match dagilim_min] = match_kDtree(q,nokta,kdOBJ);
    end

    if arg.EdgeRejection
        n_indeks = not(ismember(match, bdr));
        veri_indeks = match(n_indeks);
        dagilim_min = dagilim_min(n_indeks);
    else
        n_indeks = true(1, Np);
        veri_indeks = match;
    end

    if arg.WorstRejection
        kenar = round((1-arg.WorstRejection)*sum(n_indeks));
        ciftler = find(n_indeks);
        [~, idx] = sort(dagilim_min);
        n_indeks(ciftler(idx(kenar:end))) = false;
        veri_indeks = match(n_indeks);
        dagilim_min = dagilim_min(n_indeks);
    end
    
    if k == 1
        Hata(k) = sqrt(sum(dagilim_min.^2)/length(dagilim_min));
    end
    
    switch arg.Minimize
        case 'point'
            % agirlik vektoru
            agirlik = arg.Weight(match);
            [R,T] = eq_point(q(:,veri_indeks),nokta(:,n_indeks), agirlik(n_indeks));
        case 'plane'
            agirlik = arg.Weight(match);
            [R,T] = eq_plane(q(:,veri_indeks),nokta(:,n_indeks),arg.Normals(:,veri_indeks),agirlik(n_indeks));
        case 'lmaPoint'
            [R,T] = eq_lmaPoint(q(:,veri_indeks),nokta(:,n_indeks));
    end

    % donusume ekle
    TR(:,:,k+1) = R*TR(:,:,k);
    TT(:,:,k+1) = R*TT(:,:,k)+T;

    % son donusumu uygula
    nokta = TR(:,:,k+1) * p + repmat(TT(:,:,k+1), 1, Np);
    
    % ortalama
    Hata(k+1) = rms_error(q(:,veri_indeks), nokta(:,n_indeks));

    if arg.Extrapolation
        qq(:,k+1) = [rmat2quat(TR(:,:,k+1));TT(:,:,k+1)];
        dq(:,k+1) = qq(:,k+1) - qq(:,k);
        theta(k+1) = (180/pi)*acos(dot(dq(:,k),dq(:,k+1))/(norm(dq(:,k))*norm(dq(:,k+1))));
        if arg.Verbose
            disp(['Direction change ' num2str(theta(k+1)) ' degree in iteration ' num2str(k)]);
        end
        if k>2 && theta(k+1) < 10 && theta(k) < 10
            d = [Hata(k+1), Hata(k), Hata(k-1)];
            v = [0, -norm(dq(:,k+1)), -norm(dq(:,k))-norm(dq(:,k+1))];
            vmax = 25 * norm(dq(:,k+1));
            dv = ekstrapolasyon(v,d,vmax);
            if dv ~= 0
                q_mark = qq(:,k+1) + dv * dq(:,k+1)/norm(dq(:,k+1));
                q_mark(1:4) = q_mark(1:4)/norm(q_mark(1:4));
                qq(:,k+1) = q_mark;
                TR(:,:,k+1) = quat2rmat(qq(1:4,k+1));
                TT(:,:,k+1) = qq(5:7,k+1);
                % toplam donusum uyarlama
                nokta = TR(:,:,k+1) * p + repmat(TT(:,:,k+1), 1, Np);
                switch arg.Matching
                    case 'bruteForce'
                        [~, dagilim_min] = match_bruteForce(q,nokta);
                    case 'Delaunay'
                        [~, dagilim_min] = match_Delaunay(q,nokta,DT);
                    case 'kDtree'
                        [~, dagilim_min] = match_kDtree(q,nokta,kdOBJ);
                end
                Hata(k+1) = sqrt(sum(dagilim_min.^2)/length(dagilim_min));
            end
        end
    end
    t_vektor(k+1) = toc;
end

if not(arg.ReturnAll)
    TR = TR(:,:,end);
    TT = TT(:,:,end);
end

function [match dagilim_min] = match_bruteForce(veri, nokta)
    k = size(nokta,2);
    l = size(veri,2);    
    match = zeros(1,k);
    dagilim_min = zeros(1,k);
    for ki=1:k
        d=zeros(1,l);
        for ti=1:3
            d=d+(veri(ti,:)-nokta(ti,ki)).^2;
        end
        [dagilim_min(ki),match(ki)]=min(d);
    end
    
    dagilim_min = sqrt(dagilim_min);


function [match dagilim_min] = match_Delaunay(veri, nokta, DT)
	match = transpose(nearestNeighbor(DT, transpose(nokta)));
	dagilim_min = sqrt(sum((nokta-veri(:,match)).^2,1));

function [match dagilim_min] = match_kDtree(~, nokta, kdOBJ)
	[match dagilim_min] = knnsearch(kdOBJ,transpose(nokta));
    match = transpose(match);

function [R,T] = eq_point(veri,nokta,agirlik)

m = size(nokta,2);
n = size(veri,2);

% agirlik normallestirme
agirlik = agirlik ./ sum(agirlik);

% veri agirlik merkezi bul ve sapmalar? hesapla
veri_cizgi = veri * transpose(agirlik);
veri_isaret = veri - repmat(veri_cizgi, 1, n);
nokta_cizgi = nokta * transpose(agirlik);
nokta_isaret = nokta - repmat(nokta_cizgi, 1, m);
% agirliklari uygula
veri_isaret = veri_isaret .* repmat(agirlik, 3, 1);
N = nokta_isaret*transpose(veri_isaret); % noktalari al
% SVD
[U,~,V] = svd(N); 
R = V*diag([1 1 det(U*V')])*transpose(U);
T = veri_cizgi - R*nokta_cizgi;

function [R,T] = eq_plane(veri,nokta,i,agirlik)
i = i .* repmat(agirlik,3,1);
c = cross(nokta,i);
cn = vertcat(c,i);
C = cn*transpose(cn);
b = - [sum(sum((nokta-veri).*repmat(cn(1,:),3,1).*i));
       sum(sum((nokta-veri).*repmat(cn(2,:),3,1).*i));
       sum(sum((nokta-veri).*repmat(cn(3,:),3,1).*i));
       sum(sum((nokta-veri).*repmat(cn(4,:),3,1).*i));
       sum(sum((nokta-veri).*repmat(cn(5,:),3,1).*i));
       sum(sum((nokta-veri).*repmat(cn(6,:),3,1).*i))];
   
X = C\b;
cx = cos(X(1)); cy = cos(X(2)); cz = cos(X(3)); 
sx = sin(X(1)); sy = sin(X(2)); sz = sin(X(3)); 
R = [cy*cz cz*sx*sy-cx*sz cx*cz*sy+sx*sz;
     cy*sz cx*cz+sx*sy*sz cx*sy*sz-cz*sx;
     -sy cy*sx cx*cy];   
T = X(4:6);

function [R,T] = eq_lmaPoint(q,p)

Rx = @(a)[1     0       0;
          0     cos(a)  -sin(a);
          0     sin(a)  cos(a)];
      
Ry = @(b)[cos(b)    0   sin(b);
          0         1   0;
          -sin(b)   0   cos(b)];
      
Rz = @(g)[cos(g)    -sin(g) 0;
          sin(g)    cos(g)  0;
          0         0       1];

Rot = @(x)Rx(x(1))*Ry(x(2))*Rz(x(3));
myfun = @(x,xdata)Rot(x(1:3))*xdata+repmat(x(4:6),1,length(xdata));
options = optimset('Algorithm', 'levenberg-marquardt');
x = lsqcurvefit(myfun, zeros(6,1), nokta, veri, [], [], options);
R = Rot(x(1:3));
T = x(4:6);

function [dv] = ekstrapolasyon(v,d,vmax)

p1 = polyfit(v,d,1); % do?rusal
p2 = polyfit(v,d,2); % parabol
v1 = -p1(2)/p1(1); % 0 kesisme
v2 = -p2(2)/(2*p2(1)); % Tepe nokta polinom

if issorted([0 v2 v1 vmax]) || issorted([0 v2 vmax v1])
    disp('Parabol guncelleme');
    dv = v2;
elseif issorted([0 v1 v2 vmax]) || issorted([0 v1 vmax v2])...
        || (v2 < 0 && issorted([0 v1 vmax]))
    disp('Cizgi guncelleme');
    dv = v1;
elseif v1 > vmax && v2 > vmax
    disp('Maksimum guncelleme');
    dv = vmax;
else
    disp('Ekstra polasyon yok');
    dv = 0;
end

function ER = rms_error(p1,p2)
dsq = sum(power(p1 - p2, 2),1);
ER = sqrt(mean(dsq));

%Matris quaternion cevrimi -wikipedi

function quaternion = rmat2quat(R)

Qxx = R(1,1,:);
Qxy = R(1,2,:);
Qxz = R(1,3,:);
Qyx = R(2,1,:);
Qyy = R(2,2,:);
Qyz = R(2,3,:);
Qzx = R(3,1,:);
Qzy = R(3,2,:);
Qzz = R(3,3,:);

w = 0.5 * sqrt(1+Qxx+Qyy+Qzz);
x = 0.5 * sign(Qzy-Qyz) .* sqrt(1+Qxx-Qyy-Qzz);
y = 0.5 * sign(Qxz-Qzx) .* sqrt(1-Qxx+Qyy-Qzz);
z = 0.5 * sign(Qyx-Qxy) .* sqrt(1-Qxx-Qyy+Qzz);

quaternion = reshape([w;x;y;z],4,[]);

%Rotasyon matris quaternion cevrimi -wikipedi

function R = quat2rmat(quaternion)
q0(1,1,:) = quaternion(1,:);
qx(1,1,:) = quaternion(2,:);
qy(1,1,:) = quaternion(3,:);
qz(1,1,:) = quaternion(4,:);

R = [q0.^2+qx.^2-qy.^2-qz.^2 2*qx.*qy-2*q0.*qz 2*qx.*qz+2*q0.*qy;
     2*qx.*qy+2*q0.*qz q0.^2-qx.^2+qy.^2-qz.^2 2*qy.*qz-2*q0.*qx;
     2*qx.*qz-2*q0.*qy 2*qy.*qz+2*q0.*qx q0.^2-qx.^2-qy.^2+qz.^2];
 
function n = lsqnormest(p, k)
m = size(p,2);
n = zeros(3,m);

v = ver('stats');
if str2double(v.Version) >= 7.5 
    neighbors = transpose(knnsearch(transpose(p), transpose(p), 'k', k+1));
else
    neighbors = k_nearest_neighbors(p, p, k+1);
end

for i = 1:m
    x = p(:,neighbors(2:end, i));
    p_bar = 1/k * sum(x,2);
    
    P = (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k)); %spd matrix P
    %P = 2*cov(x);
    
    [V,D] = eig(P);
    
    [~, idx] = min(diag(D)); % choses the smallest eigenvalue
    
    n(:,i) = V(:,idx);   % returns the corresponding eigenvector    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% k-d agaci en yakin noktalar bulma
% Oklid mesafe

function [komsular mesafe] = k_nearest_neighbors(veriMatris, siraMatris, i)

veriNoktalari_sayisi = size(veriMatris,2);
siraNoktalari_sayisi = size(siraMatris,2);

komsular = zeros(i,siraNoktalari_sayisi);
mesafe = zeros(i,siraNoktalari_sayisi);

D = size(veriMatris, 1);

for i=1:siraNoktalari_sayisi
    d=zeros(1,veriNoktalari_sayisi);
    for t=1:D 
        d=d+(veriMatris(t,:)-siraMatris(t,i)).^2;
    end
    for j=1:i
        [s,t] = min(d);
        komsular(j,i)=t;
        mesafe(j,i)=sqrt(s);
        d(t) = NaN; 
    end
end

function sinir = find_bound(nokta, polinom)


polinom = double(polinom);
nokta = double(nokta);

TR = TriRep(polinom, nokta(1,:)', nokta(2,:)', nokta(3,:)');
FF = freeBoundary(TR);
sinir = FF(:,1);
 


 