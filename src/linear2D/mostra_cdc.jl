function mostra_cdc(SEGMENTOS,PONTOS,NOS,CDC,ELEM,CCSeg)
% Programa para visualização das condições de contorno em problemas
% de temperatura
%


nsegs=length(SEGMENTOS(:,1));

maxT=0;
maxq=0;

% Define os valores máximos da temperatura e do fluxo
for i=1:nsegs
    tipoCDC=CCSeg(i,2); % Tipo da condição de contorno
    valorCDC=CCSeg(i,3); % Valor da condição de contorno
    if(tipoCDC==0)
        if(maxT<abs(valorCDC))
            maxT=valorCDC;
        end
    else
        if(maxq<abs(valorCDC))
            maxq=valorCDC;
            if(valorCDC<0)
                sinalq=-1;
            else
                sinalq=1;
            end
        end
    end
end


% Definição da maior coordenada x ou y (maxp)
if max(PONTOS(:,2)) > max(PONTOS(:,3))
    maxp = max(PONTOS(:,2));
else
    maxp = max(PONTOS(:,3));
end;
fat = .01*maxp;		% Fator utilizado na construção da linha
% de divisão entre dois elementos

qsi=[-1,0,1];

for el = 1 : length(ELEM(:,1));	% for over the elements
    if (maxT ~= 0)
        fatT = CDC(el,3)/maxT;
    else fatT = 0;
    end;
    if (maxq ~= 0)
        fatq = sinalq*CDC(el,3)/maxq;
    else
        fatq = 0;
    end;
    inos = [ELEM(el,2),ELEM(el,3)];
    x = [NOS(inos(1),2),NOS(inos(2),2)];
    y = [NOS(inos(1),3),NOS(inos(2),3)];
    no1 = ELEM(el,2);
    no2 = ELEM(el,3);
    x1 = NOS(no1,2);	y1 = NOS(no1,3);
    x2 = NOS(no2,2);	y2 = NOS(no2,3);
    L=sqrt((x2-x1)^2+(y2-y1)^2); % Comprimento do elemento
    
    sx=(x2-x1)/L; % Componente x do vetor tangente
    sy=(y2-y1)/L; % Componente y do vetor tangente
    
    for noloc = 1 : 2					% for over the local nodes
        x_no=x(noloc);
        y_no=y(noloc);
        if (CDC(el,2*noloc)== 0)		% a temperatura é conhecida
            xcc = x_no+7*fat*sy*fatT;
            ycc = y_no-7*fat*sx*fatT;
            line([x_no xcc],[y_no ycc],'color',[.8 .4 .2]);
            plot(xcc,ycc,'o','color',[.6 .4 .2]);
        elseif fatq > 0
            xcc = x_no+5.5*fat*sy*fatq;
            ycc = y_no-5.5*fat*sx*fatq;
            xcc1 = xcc+1.5*fat*sy;
            ycc1 = ycc-1.5*fat*sx;
            xcc2 = xcc-fat*sx;
            ycc2 = ycc-fat*sy;
            xcc3 = xcc+fat*sx;
            ycc3 = ycc+fat*sy;
            line([x_no xcc],[y_no ycc],'color','r');
            line([xcc1 xcc2 xcc3 xcc1],[ycc1 ycc2 ycc3 ycc1],'color','r');
        elseif fatq <= 0
            xcc = x_no+7*fat*sy*abs(fatq);
            ycc = y_no-7*fat*sx*abs(fatq);
            xcc1 = x_no+1.5*fat*sy;
            ycc1 = y_no-1.5*fat*sx;
            xcc2 = xcc1-fat*sx;
            ycc2 = ycc1-fat*sy;
            xcc3 = xcc1+fat*sx;
            ycc3 = ycc1+fat*sy;
            if fatq < 0
                line([xcc1 xcc],[ycc1 ycc],'color','r');
                line([xcc2 xcc3 x_no xcc2],[ycc2 ycc3 y_no ycc2],'color','r');
            else
                line([xcc2 xcc3 x_no xcc2],[ycc2 ycc3 y_no ycc2],'color','black');
            end;
        end;
        
    end;
end;

drawnow		% Faz os gráficos pendentes

