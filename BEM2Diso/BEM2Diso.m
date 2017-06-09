%% Programa que calcula as propriedades geom�tricas de figuras planas
% usando elementos de contorno quadr�ticos cont�nuos e descont�nuos
%%
tic
addpath nrb

clear; % Apaga as vari�veis da mem�ria
close all; % Fecha todas as figuras
clc; % Limpa as linhas de comando
dad_1
%crv=dad2nurbs2(PONTOS,SEGMENTOS); % Sem parte gr�fica - Comentar linhas 87
% a 129 e 138.
[crv,contorno]=dad2nurbs(PONTOS,SEGMENTOS); % Com parte gr�fica - Retirar coment�rio
% linhas 87 a 129 e 138.
n=size(crv,2);
for i=1:n
    dcrv(i)=nrbderiv(crv(i));
end
p=1;%refinamento p (aumenta a ordem das NURBS)
for i=1:n
crv(i) = nrbdegelev(crv(i), p);    
end

h=10;%refinamento h (aumenta o n�mero de NURBS)
for i=1:n
    novosnos=linspace(0,1,h+2);
crv(i) = nrbkntins(crv(i), novosnos(2:end-1));    
end

% Cria os pontos de coloca��o
z=0;
for k=1:n
    p=crv(k).order-1;
    nnos(k)=crv(k).number;
    valorCDC=CCSeg(k,3);
    tipoCDC=CCSeg(k,2);
    for i=1:crv(k).number
        z=z+1;
        numcurva(z)=k;
        collocPts(z)=sum(crv(k).knots((i+1):(i+p)))/p;
        if(i==2)
            collocPts(z-1)=(collocPts(z)+collocPts(z-1))/2;
        end
        if(i==nnos(k))
            collocPts(z)=(collocPts(z)+collocPts(z-1))/2;
        end
        
       CDC(z,:) = [z,tipoCDC,valorCDC];       
    end
end


% Cria a matriz E para transfer�ncia das condi��es de contorno dos pontos
% de controle para os pontos de coloca��o

nnos=[0 nnos];
nnos=cumsum(nnos);

E=zeros(length(collocPts));
for i=1:length(collocPts)
    collocCoord(i,:)=nrbeval(crv(numcurva(i)), collocPts(i));
    [B, id] = nrbbasisfun (collocPts(i),crv(numcurva(i)));
    E(i,id+1+nnos(numcurva(i)))=B;
end

% Plota os pontos de coloca��o
figure
plot(collocCoord(:,1),collocCoord(:,2),'rd')
for i=1:n
    nrbplot(crv(i),50)%plota a figura
    hold on
    plot(crv(i).coefs(1,:)./crv(i).coefs(4,:),crv(i).coefs(2,:)./crv(i).coefs(4,:),'r--')%plota o pol�gono de controle
    hold on
    plot(crv(i).coefs(1,:)./crv(i).coefs(4,:),crv(i).coefs(2,:)./crv(i).coefs(4,:),'*k')%plota os pontos de controle
    hold on
end
pg=12;% N�mero de pontos de Gauss usado na integra��o
% Cria as matriz as matrizes H e G
[H,G]=calcula_PropGeom(n,collocCoord,collocPts,nnos,crv,dcrv,pg,kmat); % Calcula as propriedades geom�tricas de figuras planas
H=H+E/2;

% Aplica as condi��es de contorno
[A,b]= aplica_CDC(G,H,CDC,E);
x=A\b; % Calcula o vetor x

[Tc,qc]=monta_Teq(CDC,x); % Separa temperatura e fluxo

%% Gera��o da malha
% C�lculo das integrais nas trimmed surfaces
ielem=0;
subd=17;
u=linspace(0.0,1.0,subd);
ncont=size(contorno,1);
ELEM=[1 2];
for j=1:ncont
    pini=contorno(j,1); % Ponto onde come�a o contorno j
    ncurv=contorno(j,2);
    ielemp=ielem;
    for i=1:ncurv
        p=nrbeval(crv(pini+i-1),u);
        %Range
        numP=size(p,2);
        NOS(1+ielem:numP+ielem,1:3)=p';
        ELEM(1+ielem:numP+ielem,1)=(1+ielem:numP+ielem)';
        ELEM(1+ielem:numP+ielem,2)=(2+ielem:(numP+1+ielem))';
        ielem=ielem+numP;
    end
    ELEM(ielem,2)=ielemp+1;
end

nelem=size(ELEM,1);
xymax=max(NOS);
xymin=min(NOS);

% Triangula��o de Delaunay
% t=tri�ngulos

t = delaunay(NOS(:,1:2)); 
%t = delaunay(collocCoord(:,1:2)); 

% Impose geometry constraints
% Triangle centroids
centroid=(NOS(t(:,1),1:2)+NOS(t(:,2),1:2)+NOS(t(:,3),1:2))/3.0;
% Verifica se os centroides est�o dentro ou fora do dom�nio
i = inpoly(centroid,NOS(:,1:2),ELEM);                                      
% Take triangles with internal centroids
t = t(i,:); % Somente tri�ngulos internos

% Plota os tri�ngulos para mostrar a trimmed surface
patch('faces',t,'vertices',NOS,'FaceColor','b');
%%

% Transfere temperatura e fluxo dos pontos de controle para os pontos de
% coloca��o
T=E*Tc';
q=E*qc';

%% Mapa de cores
mapa_de_cor(NOS,t,T)