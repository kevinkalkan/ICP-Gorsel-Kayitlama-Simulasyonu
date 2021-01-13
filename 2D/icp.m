function [Rotasyon_T,Cevrilmis_T,veri] = icp(model,veri,Max_I,Min_I,krit,esik)

%   model - olusturulmus olan son nokta bulutu
%   veri - rastgele alinan verilerin matrisi
%   R - Rotasyon matrisi
%   T - Cevrim vektoru
%   data2 - donusturulmus veri noktalari matrisi

if nargin <6
    esik=1e-5;                     % iterasyon esigi
    if nargin<5
        krit=0;                  
        if nargin<4
            Min_I=15;              %iterasyon 
            if nargin<3
                Max_I=100;        
            end
        end
    end
    
end

if or(isempty(model),isempty(veri))
    error('Model ve veri noktalarinda problem');
end

if isempty(Max_I)
    Max_I=200;
end

if isempty(Min_I)
    Min_I=10;
end

if isempty(krit)
    krit=0;
end

if isempty(esik)
    esik=1e-5;
end

% Noktalarin boyutlari

if (size(model,2)<size(model,1))
    Transpoz_matris=true;
    matris=size(model,2);
    Matris=size(model,1);
else
    Transpoz_matris=false;
    matris=size(model,1);
    Matris=size(model,2);
end
if (size(veri,2)<size(veri,1))
    veri=veri';
end
N=size(veri,2);

% En yakin nokta arama yapisi 

if matris<4
    if Transpoz_matris
        DT=delaunayTriangulation(model);
    else
        DT=delaunayTriangulation(model');
    end
else
    DT=[];
    resid=zeros(N,1);
    inp=ones(N,1);
end

% Dönüstürme

Rotasyon_T=eye(matris);
Cevrilmis_T=zeros(matris,1);

%ICP

ayar=9e99;

for iter=1:Max_I   
    eski_ayar=ayar;  
    % Veri noktalarina en yakin son model noktalari bulma  
    if isempty(DT)
        if Transpoz_matris
            for i=1:N
                min_deger=9e99;
                for j=1:Matris
                    val=norm(veri(:,i)-model(j,:)');
                    if val<min_deger
                        min_deger=val;
                        inp(i)=j;
                        resid(i)=val;
                    end
                end
            end
        else
            for i=1:N
                min_deger=9e99;
                for j=1:Matris
                    val=norm(veri(:,i)-model(:,j));
                    if val<min_deger
                        min_deger=val;
                        inp(i)=j;
                        resid(i)=val;
                    end
                end
            end
        end
    else
        [inp,resid] = nearestNeighbor(DT,veri');
    end  
    % Donusum bulma
    switch krit
            case 0           
            ayar=mean(resid.^2);
            med=mean(veri,2);
            if Transpoz_matris
                mem=mean(model(inp,:),1);
                C=veri*model(inp,:)-(N*med)*mem;
                [U,~,V]=svd(C);
                Ri=V*U';
                if det(Ri)<0
                    V(:,end)=-V(:,end);
                    Ri=V*U';
                end
                Ti=mem'-Ri*med;
            else
                mem=mean(model(:,inp),2);
                C=veri*model(:,inp)'-(N*med)*mem';
                [U,~,V]=svd(C);
                Ri=V*U';
                if det(Ri)<0
                    V(:,end)=-V(:,end);
                    Ri=V*U';
                end
                Ti=mem-Ri*med;
            end
    end
    
    veri=Ri*veri;  % Donusum yap
    for i=1:matris
        veri(i,:)=veri(i,:)+Ti(i);      
    end
    Rotasyon_T=Ri*Rotasyon_T;    % Donusum guncelleme
    Cevrilmis_T=Ri*Cevrilmis_T+Ti;                        
    if iter >= Min_I
        if abs(eski_ayar-ayar) < esik
            break
        end
    end
    
end

