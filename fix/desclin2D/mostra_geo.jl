function mostra_geo(SEGMENTOS,PONTOS,ELEM,NOS)

% Programa de visualização da geometria e malha
nelem=length(ELEM(:,1));
figure
axes('DataAspectRatio',[1,1,1]);
% axis off;
hold on

% Definition of the biggest dimension of the problem
max_dl = 0;
for lin = 1 : length(SEGMENTOS(:,1))
    p1 = SEGMENTOS(lin,2);
    p2 = SEGMENTOS(lin,3);
    xp1 = PONTOS(p1,2);	yp1 = PONTOS(p1,3);
    xp2 = PONTOS(p2,2);	yp2 = PONTOS(p2,3);
    dl = sqrt((xp1-xp2)^2+(yp1-yp2)^2);
    if dl > max_dl
        max_dl = dl;
    end;
end;


fat = .02*max_dl;		% Scale factor

% Plot two PONTOS in white color only to define the size of the dimension 
% of the window in which the geometry will be shown

emin = min(PONTOS(:,2))-.1*max_dl;
dmax = max(PONTOS(:,2))+.1*max_dl;
imin = min(PONTOS(:,3))-.1*max_dl;
smax = max(PONTOS(:,3))+.1*max_dl;
plot(emin,imin,'w');	% Botton left point
plot(dmax,smax,'w');	% Top right point

% Plot the PONTOS that defines the geometry as a black x.
for i = 1 : length(PONTOS(:,1))
    plot(PONTOS(i,2),PONTOS(i,3),'blackx','markersize',8);
end;

% Generate PONTOS over elements that will be used to interpolate geometry
ponti = -1 : .2 : 1;
qsi_ponta=[-1,1];
% Compute the continuous shape functions
for el = 1 : nelem
    
    
    no1 = ELEM(el,2);
    no2 = ELEM(el,3);
    
    for p = 1 : 11
        [N1,N2] = calc_fforma(ponti(p));
        FF(p,1)=N1;
        FF(p,2)=N2;
    end;
    
    
    
    x1 = NOS(no1,2);	y1 = NOS(no1,3);
    x2 = NOS(no2,2);	y2 = NOS(no2,3);
    L=sqrt((x2-x1)^2+(y2-y1)^2); % Comprimento do elemento
    
    sx=(x2-x1)/L; % Componente x do vetor tangente
    sy=(y2-y1)/L; % Componente y do vetor tangente
    n1=sy; % Componente x do vetor normal
    n2=-sx; % Componente y do vetor normal
    
    % Compute the PONTOS to interpolate the geometry over the element
    x_el = [x1;x2];
    y_el = [y1;y2];
    XC = FF*x_el;	% Vector (10x1) with the x coordinates of the PONTOS
    YC = FF*y_el;	% Vector (10x1) with the y coordinates of the PONTOS
    
    
    % Plot the element
    plot(XC,YC,'color','black','LineWidth',1.2,'LineStyle','-')
    plot(x_el,y_el,'k.','markersize',6);	% Plot the node of the elements
    % Plot a line in the beginning and in the end of the element
    plot([XC(1)+fat*n1 XC(1)-fat*n1], ...
        [YC(1)+fat*n2(1) YC(1)-fat*n2],'LineWidth',1.2,'color', ...
        'black','LineStyle','-');
    plot([XC(end)+fat*n1 XC(end)-fat*n1], ...
        [YC(end)+fat*n2 YC(end)-fat*n2],'LineWidth',1.2,'color',...
        'black','LineStyle','-');
end;

