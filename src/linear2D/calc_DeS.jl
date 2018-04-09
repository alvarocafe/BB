function [Dx,Dy,Sx,Sy]=calc_DeS(x1,y1,x2,y2,xd,yd,k,qsi,w)

npg=length(qsi); % Número de pontos de Gauss

Dx=[0 0]; % Inicializa a matriz h do elemento
Dy=[0 0]; % Inicializa a matriz g do elemento
Sx=[0 0]; % Inicializa a matriz h do elemento
Sy=[0 0]; % Inicializa a matriz g do elemento

L=sqrt((x2-x1)^2+(y2-y1)^2); % Comprimento do elemento
dgamadqsi=L/2; % Jacobiano

sx=(x2-x1)/L; % Componente x do vetor tangente
sy=(y2-y1)/L; % Componente y do vetor tangente
nx=sy; % Componente x do vetor normal
ny=-sx; % Componente y do vetor normal

for kk=1:npg
    [N1,N2]=calc_fforma(qsi(kk)); % Calcula as funções de forma

    x=N1*x1+N2*x2; % Calcula a coordenada x do ponto de integração
    y=N1*y1+N2*y2; % Calcula a coordenada y do ponto de integração
    
    [dTdx,dTdy,dqdx,dqdy]=calc_dsolfund(x,y,xd,yd,nx,ny,k); % Calcula 
       % as derivadas das soluções fundamentais
    Dx=Dx+dTdx*[N1,N2]*dgamadqsi*w(kk); % Integral da matriz H
    Dy=Dy+dTdy*[N1,N2]*dgamadqsi*w(kk); % Integral da matriz H
    Sx=Sx+dqdx*[N1,N2]*dgamadqsi*w(kk); % Integral da matriz G
    Sy=Sy+dqdy*[N1,N2]*dgamadqsi*w(kk); % Integral da matriz G
   
end