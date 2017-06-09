function crv=dad2nurbs2(PONTOS,SEGMENTOS)

contorno=calc_ncont(SEGMENTOS);
ncont=size(contorno,1);
icrv=0;
for j=1:ncont
%     icolA=1;
    pini=contorno(j,1); % Ponto onde começa o contorno j
    nseg=contorno(j,2);
%     temcirc=0;
    for i=1:nseg
        icrv=icrv+1;
        raio=SEGMENTOS(i+pini-1,4);%define valores para o raio
        np1=SEGMENTOS(i+pini-1,2); % Define as coordenas x
        np2=SEGMENTOS(i+pini-1,3); %Define as coordenasdas y
        p1=[PONTOS(np1,2), PONTOS(np1,3),0]; %Define o primeiro ponto da curva
        p2=[PONTOS(np2,2), PONTOS(np2,3),0]; %Define o segundo ponto da curva
        if(raio==0)
            crv(icrv)=nrbline(p1,p2);
        else
            [xc,yc]=calcula_centro(p1(1),p1(2),p2(1),p2(2),raio);
            [sang,eang] = calcula_arco(p1(1),p1(2),p2(1),p2(2),xc,yc,raio);
%             sang=atan2(p1(2)-yc,p1(1)-xc);
%             eang=atan2(p2(2)-yc,p2(1)-xc);
            crv(icrv) = nrbcirc(raio,[xc,yc,0],sang,eang);
%             ncoefs=size(crv(i).coefs,2);%número de linhas de coefs
           
%             temcirc=1;
%              if SEGMENTOS(:,4)==zeros(4,1)
%             crv(i)=cr;
%             
%             crv
        end
 % nrbplot(crv(i),50);
 % hold on;   
    end

 
end

return