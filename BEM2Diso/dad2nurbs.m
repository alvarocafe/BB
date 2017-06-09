function [crv,contorno]=dad2nurbs(PONTOS,SEGMENTOS)

contorno=calc_ncont(SEGMENTOS);
ncont=size(contorno,1);
icrv=0;
jj=0;
for j=1:ncont   
    pini=contorno(j,1); % Ponto onde começa o contorno j
    nseg=contorno(j,2);
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
            
            crv(icrv) = nrbcirc(raio,[xc,yc,0],sang,eang);
        end
%        jj=jj+1;
%        sup(jj) = nrbextrude(crv(icrv),[0 0 espessura]);        
    end  
end
return