function format_dad(PONTOS,SEGMENTOS,MALHA,CCSeg)
    #Programa para formatar os dados de entrada
    #Autor: Ãlvaro Campos Ferreira, Ã‰der Lima de Albuquerque

    #Vamos declarar as variÃ¡veis que estÃ£o com type instabilities
    r1::Float64 = 0;
    r2::Float64 = 0;
    tet::Float64 = 0;
    tet2::Float64 = 0;
    sig::Int64 = 1;
    divtet::Float64 = 0;
    x_i::Float64 = 0;
    y_i::Float64 = 0;
    x_m::Float64 = 0;
    y_m::Float64 = 0;
    x_f::Float64 = 0;
    y_f::Float64 = 0;

    #Definimos o tamanho das matrizes
    tamanho::Int64 = 0   #Iniciamos o tamanho

    for i = 1 : size(MALHA)[1]
        tamanho = tamanho + MALHA[i,2]
    end
    NOS_GEO = zeros(tamanho,3)
    ELEM = zeros(Int,tamanho,3)
    NOS = zeros(tamanho,3)
    normal = zeros(tamanho,tamanho)
    cont_nos = 0;
    cont_el = 0;
    num_lin = size(SEGMENTOS)[1]
    p_ini::Int64 = SEGMENTOS[1,2]
    #DefiniÃ§Ã£o da maior dimensÃ£o do problema
    max_dl::Float64 = 0
    for lin = 1 : num_lin
        p1::Int64 = SEGMENTOS[lin,2]
        p2::Int64 = SEGMENTOS[lin,3]
        xp1 = PONTOS[p1,2]
        yp1 = PONTOS[p1,3]
        xp2 = PONTOS[p2,2]
        yp2 = PONTOS[p2,3]
        dl = sqrt((xp1 - xp2)^2 + (yp1 - yp2)^2)
        if dl > max_dl
            max_dl = dl
        end
    end

    no_ini = 1
    t = 1
    pp2::Int64 = 2
    no1_prox = 0
    while(t<num_lin)  	# While over all lines
        while(pp2!=p_ini)
            num_el_lin = MALHA[t,2];	# Number of the elements in the line t
            # Coordinates of the initial and final PONTOS of each line
            # (x1l,y1l,x2l,y2l)
            pp1::Int64  = SEGMENTOS[t,2];
            pp2 = SEGMENTOS[t,3];
            x1l = PONTOS[pp1,2];
            y1l = PONTOS[pp1,3];
            x2l = PONTOS[pp2,2];
            y2l = PONTOS[pp2,3];
            # 1. Generation of the matrices NOS, NOS_GEO, NOS_DRM, ELEM e ELEM_GEO
            if(SEGMENTOS[t,4]==0) # The segment is a straight line
                # Increment in x and y direction
                delta_x = x2l - x1l;
                delta_y = y2l - y1l;
            else #The segment is an arc
                # Compute the center of the arc and its coordinates
                r = SEGMENTOS[t,4];
                xc, yc = calcula_centro(x1l,y1l,x2l,y2l,r);
                # Distance between p1 and c (r1) and between p2 and c (r2)
                r1 = sqrt((x1l-xc)^2+(y1l-yc)^2);
                r2 = sqrt((x2l-xc)^2+(y2l-yc)^2);
                if abs(r1-r2)<.00001*max_dl
                    # Compute the angle between the lines from point c to p1 (tet1) and c to p2 (tet2)
                    tet1, tet2 = calcula_arco(x1l,y1l,x2l,y2l,xc,yc);
                    if tet2 < tet1
                        tet2 = tet2 + 2*pi;
                    end;

                    # Angle of the sector defined by the arc
                    if SEGMENTOS[t,4] > 0
                        tet = abs(tet2-tet1);
                        sig = 1;
                    else
                        tet = 2*pi-abs(tet2-tet1);
                        sig = -1;
                    end;

                    # Angle between two nodes of the line
                    divtet = tet/(2*num_el_lin);
                else
                    println("Error in the data input file: Wrong central point");
                end;
            end
            # Generation of elements and nodes
            for i = 1 : num_el_lin

                if(SEGMENTOS[t,4]==0) # The segment is a straight line
                    x_i = x1l + delta_x/num_el_lin*(i-1);			# initial x coordinate of the element
                    y_i = y1l + delta_y/num_el_lin*(i-1);			# initial y coordinate of the element
                    x_m = x1l + delta_x/num_el_lin*(i-.5);	# midpoint x coordinate of the element
                    y_m = y1l + delta_y/num_el_lin*(i-.5);	# midpoint y coordinate of the element
                    lx=x_m-x_i;                              # distance in x direction between geometric nodes 1 and 2
                    ly=y_m-y_i;                              # distance in y direction between geometirc nodes 1 and 2

                else  # The segment is an arc
                    # Compute the node coordinates
                    x_i = xc+r1*cos(tet1+2*(i-1)*sig*divtet);
                    y_i = yc+r1*sin(tet1+2*(i-1)*sig*divtet);
                    x_m = xc+r1*cos(tet1+(2*i-1)*sig*divtet);
                    y_m = yc+r1*sin(tet1+(2*i-1)*sig*divtet);
                end;
                cont_el = cont_el + 1;
                # Set the coordinate of the physical and DRM nodes
                if(no1_prox==0) # the first node needs to be created
                    cont_nos = cont_nos + 1;
                    NOS_GEO[cont_nos,:]=[cont_nos,x_i,y_i];
                    no1=cont_nos;
                else
                    no1=no1_prox;
                end;
                if(pp2!=p_ini || i<num_el_lin)
                    cont_nos = cont_nos + 1;
                    if(SEGMENTOS[t,4]==0) # Straight line
                        x_f=x_m+lx;
                        y_f=y_m+ly;
                    else                 # arc
                        x_f = xc+r1*cos(tet1+(2*i)*sig*divtet);
                        y_f = yc+r1*sin(tet1+(2*i)*sig*divtet);
                    end;
                    NOS_GEO[cont_nos,:]=[cont_nos,x_f,y_f];
                    no3=cont_nos;
                    no1_prox=no3;
                else
                    if(no_ini==0)
                        cont_nos = cont_nos + 1;
                        if(SEGMENTOS[t,4]==0) # Straight line
                            x_f=x_m+lx;
                            y_f=y_m+ly;
                        else                 # Arc
                            x_f = xc+r1*cos(tet1+(2*i)*sig*divtet);
                            y_f = yc+r1*sin(tet1+(2*i)*sig*divtet);
                        end;
                        NOS_GEO[cont_nos,:]=[cont_nos,x_f,y_f];
                        no3=cont_nos;
                        no1_prox=0;
                    else
                        no3=no_ini;
                        no1_prox=0;
                    end;
                end;
                ELEM[cont_el,:]=[cont_el,no1,no3];
            end;# end of for i = 1 : num_el_lin
            if pp2 == p_ini
                if t < num_lin
                    p_ini = SEGMENTOS[t+1,2];
                    if(SEGMENTOS[t+1,3]==2)
                        no_ini = 0;
                    else
                        no_ini = cont_nos+1;
                    end
                end;
            end;
t=t+1;
end;                                  #end of while p2
end

# Geraï¿½ï¿½o da matriz CDC (Condiï¿½ï¿½es de Contorno)
# CDC = [n. do elemento, tipo de cdc, valor da cdc]
# Tipos de cdc: 0 : temperatura conhecida
#               1 : fluxo conhecido
n_elem=size(ELEM,1);
cont_el2 = 0;
CDC=zeros(n_elem,5);
for l = 1 : length(SEGMENTOS[:,1])
    n_el_lin = MALHA[l,2];
    el_ini::Int64 = cont_el2 + 1;
    el_fin::Int64 = cont_el2 + n_el_lin;
    valorCDC=CCSeg[l,3];
    # valorCDC=complex(CCSeg[l,3],CCSeg[l,4]);
    tipoCDC=CCSeg[l,2];
    for el = el_ini : el_fin
        CDC[el,:] = [el,tipoCDC,valorCDC,tipoCDC,valorCDC];
    end;
    cont_el2 = el_fin;
end;

nnos=size(ELEM,1); # Nï¿½mero de nï¿½s
NOS=zeros(nnos,3);
for i=1:nnos # Laï¿½o sobre os pontos fontes
    pontoi::Int64=ELEM[i,2]; # Ponto final do elemento
    pontof::Int64=ELEM[i,3]; # Ponto inicial do elemento
    xi=NOS_GEO[pontoi,2]; # Coordenada x do ponto inicial do elemento
    xf=NOS_GEO[pontof,2]; # Coordenada x do ponto final do elemento
    yi=NOS_GEO[pontoi,3]; # Coordenada y do ponto inicial do elemento
    yf=NOS_GEO[pontof,3];  # Coordenada y do ponto final do elemento
    xd=(xi+xf)/2; # Coordenada x do ponto fonte
    yd=(yi+yf)/2; # Coordenada y do ponto fonte
    h1=xf-xi;
    h2=yf-yi;
    el = sqrt(h1^2 + h2^2);
    normal[1,i] = h2/el;
    normal[2,i] = -h1/el;
    NOS[i,1]=i;
    NOS[i,2]=xd;
    NOS[i,3]=yd;
end
return NOS_GEO,NOS,ELEM,CDC
end
function calc_fforma_d(qsi)
    # Evaluates the discontinuous linear shape functions for points located at a quarter from the endpoints of the element.
    N1=1/2 -1*qsi;
    N2=1/2 + 1*qsi;
    return N1,N2
end
function	calc_fforma(qsi)
    # Evaluates the continuous linear shape function.
    N1=1/2*(1-qsi);
    N2=1/2*(1+qsi);
  
    return N1,N2
end
% Função para testar a proximidade de um ponto interno ao contorno
% Se o ponto estiver mais proximo que o desejado (d<d_min) o ponto
% é reprovado (aceita = não). Caso contrrário (aceita = sim)
%
%   Autor: Frederico Lourenço
%   Data de criação setembro de 1999 
%   Revisão 0.0

function [aceita] = aceita_dist(xpi,ypi,xl1,yl1,xl2,yl2,d_min,raio)
  
    pc = 0;					% Contador de pontos do contorno
        npc = length(xl1);	% N. de pontos do contorno
            aceita = 'sim';
          
            % Verificando a proximidade do ponto interno aos pontos do contorno		
                while (strcmp(aceita,'sim') && pc < npc)
                    pc = pc + 1;
                    x1 = xl1(pc);
                    y1 = yl1(pc);
                    d = sqrt((xpi-x1)^2+(ypi-y1)^2);
                    if d < d_min
                        aceita = 'não';
                    end;
                end;
              
                % Verificando a proximidade às linhas do contorno
                    l = 0;	% Contador das linhas do contorno
                      
                        while (strcmp(aceita,'sim') && l < npc)
                            l = l + 1;
                            x1 = xl1(l); y1 = yl1(l);
                            x2 = xl2(l); y2 = yl2(l);   
                            if(raio(l)==0) % The segment is a straight line
                              
                                if (xpi > x1 && xpi < x2) || (xpi > x2 && xpi < x1)
                                    m = (y2-y1)/(x2-x1);		% Angular coeficient of the recta
                                    yi = m*(xpi-x1)+y1;			% y coordinate of the intersection of a recta normal to the current boundary
                                    % recta that cross the point
                                    dy = ypi - yi; 
                                    d = abs(dy*cos(atan(m))); % distance from the point to the recta
                                    if d < d_min 
                                        aceita = 'não';    
                                    end;
                                end;
                              
                                if (ypi > y1 && ypi < y2) || (ypi > y2 && ypi < y1)
                                    if x1 == x2
                                        d = abs(xpi-x1);
                                    else 
                                        m = (y2-y1)/(x2-x1);			% Angular coeficient of the recta
                                        xi = 1/m*(ypi-y1)+x1;	% x coordinate of the intersection of a recta normal to the current boundary
                                        % recta that cross the point            
                                        dx = xpi - xi;
                                        d = abs(dx*sin(atan(m))); % distance from the point to the recta
                                    end;
                                    if d < d_min
                                        aceita = 'não';    
                                    end;
                                end;
                            else    % The segment is an arc
                                [xc,yc]=calcula_centro(x1,y1,x2,y2,raio(l)); % Center of the arc
                                [teta_i,teta_f]=calcula_arco(x1,y1,x2,y2,xc,yc); % angle of the lines that defines the arc with the horizontal direction
                                [teta_i,teta_p]=calcula_arco(x1,y1,xpi,ypi,xc,yc);  % teta_p angle of the line that cross the center point and the
                                % internal point with the horizontal direction
                                if(raio(l)>0) % The center is in the left side of the arc (from the initial to the end point)
                                    if(teta_f>teta_i) 
                                        if(teta_p>teta_i && teta_p<teta_f)
                                            d=abs(raio(l)-sqrt((xpi-xc)^2+(ypi-yc)^2)); % distance from the point to the arc
                                            if d < d_min
                                                aceita = 'não';    
                                            end;
                                        end;
                                    else % teta_f<teta_i
                                        if(teta_p>teta_i || teta_p<teta_f)
                                            d=abs(raio(l)-sqrt((xpi-xc)^2+(ypi-yc)^2)); % distance from the point to the arc
                                            if d < d_min
                                                aceita = 'não';    
                                            end;
                                        end;
                                    end;
                                else % raio(l) < 0 % The center is in the right side of the arc (from the initial to the end point)
                                    if(teta_i > teta_f)
                                        if(teta_p>teta_f && teta_p<teta_i)
                                            d=abs(abs(raio(l))-sqrt((xpi-xc)^2+(ypi-yc)^2)); % distance from the point to the arc
                                            if d < d_min
                                                aceita = 'não';    
                                            end;
                                        end;
                                    else % teta_i < teta_f
                                        if(teta_p>teta_f || teta_p<teta_i)
                                            d=abs(abs(raio(l))-sqrt((xpi-xc)^2+(ypi-yc)^2)); % distance from the point to the arc
                                            if d < d_min
                                                aceita = 'não';    
                                            end;
                                        end;
                                    end;
                                end;         
                            end;
                        end;


% Retorna a aprovação ou não do ponto interno
    function calcula_arco(x1,y1,x2,y2,xc,yc)
        # Function to compute the tet1 angle between the line from point (x1,y1) to (xc,yc) and the
        # horizontal direction and the the tet2 angle between the line from point (x2,y2) to (xc,yc)
        # and the horizontal direction


        dx1 = x1 - xc; dy1 = y1 - yc;
        dx2 = x2 - xc; dy2 = y2 - yc;

        # Computation of tet1
        if dy1 == 0				# The point 1 and the center have the same y coordinate
            if x1 > xc
                tet1 = 0;
            else  # (x1 < xc)
                tet1 = pi;
            end;
        elseif dx1 == 0				# The point 1 and the center have the same x coordinate
            if y1 > yc
                tet1 = pi/2;
            else  # (y1 < yc)
                tet1 = -pi/2;
            end;
        else  # (dx1~=0 e dy1~=0)
            tet1 = atan(dy1/dx1);
            if dx1<0 && tet1<0
                tet1 = pi + tet1;
            elseif dx1 < 0 && tet1>0
                tet1 = -pi + tet1;
            end;
        end;

        # Computation of tet2
        if dy2 == 0				# The point 2 and the center have the same y coordinate
            if x2 > xc
                tet2 = 0;
            else  # (x2 < xc)
                tet2 = pi;
            end;
        elseif dx2 == 0				# The point 2 and the center have the same x coordinate
            if y2 > yc
                tet2 = pi/2;
            else  # (y2 < yc)
                tet2 = -pi/2;
            end;
        else  # (dx2~=0 e dy2~=0)
            tet2 = atan(dy2/dx2);
            if dx2<0 && tet2<0
                tet2 = pi + tet2;
            elseif dx2 < 0 && tet2>0
                tet2 = -pi + tet2;
            end;
        end;
        return tet1, tet2
    end
    function calcula_centro(x1,y1,x2,y2,raio)

        # Compute the center of an arc given two points and the radius

        xm=(x1+x2)/2;
        ym=(y1+y2)/2;
        b=sqrt((x2-x1)^2+(y2-y1)^2);
        t1=(x2-x1)/b;
        t2=(y2-y1)/b;
        n1=t2;
        n2=-t1;
        h=sqrt(abs(raio^2-(b/2)^2));
        if(raio>0)
            if(n1==0)
                xc=xm;
                yc=ym-n2/abs(n2)*h;
            else
                xc=-n1/abs(n1)*sqrt(h^2*n1^2/(n1^2+n2^2))+xm;
                yc=n2/n1*(xc-xm)+ym;
            end;
        else
            if(n1==0)
                xc=xm;
                yc=ym+n2/abs(n2)*h;
            else
                xc=n1/abs(n1)*sqrt(h^2*n1^2/(n1^2+n2^2))+xm;
                yc=n2/n1*(xc-xm)+ym;
            end;
        end;
        return xc, yc
    end
    function PONTOS_INT=gera_p_in(NPX,NPY,PONTO,SEGMENTOS)
        % Programa para criação de pontos internos a uma  geometria
        % genérica formada por retas
        %
        %   Autor: Frederico Lourenço
        %   Data de criação setembro de 1999
        %   Revisão 0.0
      
        % Definição da área máxima para criação de pontos internos
        xmin = min(PONTO(:,2));
        xmax = max(PONTO(:,2));
        ymin = min(PONTO(:,3));
        ymax = max(PONTO(:,3));
        lx = xmax - xmin;		% Largura do retângulo que contém a geometria
            ly = ymax - ymin;		% Altura do retângulo que contém a geometria
                n_SEGMENTOSs = length(SEGMENTOS(:,1));
              
                % Definição da maior SEGMENTOS do problema
                    max_dl = 0;
                    for lin = 1 : length(SEGMENTOS(:,1))
                        p1 = SEGMENTOS(lin,2);
                        p2 = SEGMENTOS(lin,3);
                        xp1 = PONTO(p1,2);	yp1 = PONTO(p1,3);
                        xp2 = PONTO(p2,2);	yp2 = PONTO(p2,3);
                        dl = sqrt((xp1-xp2)^2+(yp1-yp2)^2);
                        if dl > max_dl
                            max_dl = dl;
                        end;
                    end;
                  
                  
                    d_min = 0.03*max_dl;	% Distância mínima dos pontos internos ao contorno
                    npx = NPX+1;				% N. de pontos na horizontal
                    npy = NPY+1;				% N. de pontos na vertical
                  
                    PONTOS_INT =[];
                    % Atribuição dos pontos finais e iniciais das SEGMENTOSs aos
                    % vetores xl1, xl2, yl1 e yl2
                    for t = 1 : n_SEGMENTOSs		% Percorre todas as SEGMENTOSs
                        xl1(t) = PONTO(SEGMENTOS(t,2),2);
                        xl2(t) = PONTO(SEGMENTOS(t,3),2);
                        yl1(t) = PONTO(SEGMENTOS(t,2),3);
                        yl2(t) = PONTO(SEGMENTOS(t,3),3);
                        raio(t)= SEGMENTOS(t,4);    
                    end;
                  
                    npi = 0;	% Inicialização no n. de pontos internos
                  
                    for i = 1 : npy
                        % Criação do candidato a ponto interno (xpi,ypi)
                            ypi = ymin + (ly/npy)*i;	% y dentro do retângulo
                                for j = 1 : npx
                                    xpi = xmin + (lx/npx)*j;	% x dentro do retângulo
                                      
                                        % Início dos testes para validação do ponto interno
                                          
                                            % 1. Verificando se o ponto está dentro da geometria
                                            [xpi,ponto] = testa_ponto(xpi,ypi,xl1,yl1,xl2,yl2,lx,raio);
                                            % 2. Verificando se o ponto está muito próximo do contorno
                                                if strcmp(ponto,'interno')
                                                    aceita = aceita_dist(xpi,ypi,xl1,yl1,xl2,yl2,d_min,raio);
                                                else
                                                    aceita = 'não';
                                                end;
                                              
                                                % Armazenando os dados do ponto interno
                                                    if strcmp(aceita,'sim')	% O ponto está dentro da geometria e
                                                        npi = npi + 1;		% e respeita a distância ao contorno
                                                        PONTOS_INT(npi,:) = [npi,xpi,ypi];
                                                        plot(xpi,ypi,'.','color',[0 .5 0],'markersize',4);
                                                    end;
                                                end;
                                            end;
                                            drawnow		% Realiza as operações gráficas pendentes
                                          
                                            % Nesse ponto estão calculados os pontos internos
                                          
                                          
                                          
                                          
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
                                                                                      
                                                                                        % Mostrar numeração dos nós
                                                                                        for no  = 1 : length(NOS(:,1))
                                                                                            n_no = num2str(NOS(no,1));
                                                                                            x_no = NOS(no,2);
                                                                                            y_no = NOS(no,3);
                                                                                            text(x_no,y_no,n_no);
                                                                                        end;
                                                                                        function [xpi,ponto] = testa_ponto(xpi,ypi,xl1,yl1,xl2,yl2,lx,raio)
                                                                                            % Rotina para verificação da condição de um ponto com relação
                                                                                            % à geometria (interno ou externo).
                                                                                            % O algoritmo utilizado é baseado no número de intersecções entre
                                                                                            % o contorno e uma linha vertical que parte do ponto interno.
                                                                                                % N. ímpar de intersecções - ponto interno
                                                                                                % N. par de intersecções   - ponto externo
                                                                                                %
                                                                                                %   Autor: Frederico Lourenço
                                                                                                %   Data de criação setembro de 1999 
                                                                                                %   Revisão 0.0
                                                                                              
                                                                                              
                                                                                                pc = 0;					% Contador de pontos do contorno
                                                                                                    npc = length(xl1);	% N. de pontos do contorno
                                                                                                        sai = 'não';
		                                                                                      
                                                                                                        while (strcmp(sai,'não') && pc < npc)
                                                                                                            pc = pc + 1;
                                                                                                            if xpi == xl1(pc)
                                                                                                                xpi = xpi + lx*10^(-2);
                                                                                                                sai = 'sim';
                                                                                                            end;
                                                                                                        end;
                                                                                                      
                                                                                                        interv = 0;		% N. de intersecções entre o contorno e a linha
					                                                                % vertical abaixo do ponto interno.
                                                                                                            l=1;
                                                                                                            while l<= npc; % for over the lines that defines the boundary
                                                                                                                x1 = xl1(l); y1 = yl1(l); 
                                                                                                                x2 = xl2(l); y2 = yl2(l);
                                                                                                                if(raio(l)==0) % The segment is a straight line
                                                                                                                    if (xpi > x1 && xpi < x2) || (xpi > x2 && xpi < x1) % Check only if the x coordinate of the point is
                                                                                                                        % between the x coordinate of the intial and final points of the recta
                                                                                                                        m = (y2-y1)/(x2-x1);		% Angular coeficient of the recta
                                                                                                                        yi = m*(xpi-x1)+y1;			% y coordinate in the intersection
                                                                                                                        if ypi  >yi % If this is true, the point is above the line
                                                                                                                            interv = interv + 1;	% Counter the intersection
                                                                                                                        end;
                                                                                                                    end;
                                                                                                                    else   % The segment is an arc
                                                                                                                        [xc,yc]=calcula_centro(x1,y1,x2,y2,raio(l)); % compute the center of the arc
                                                                                                                        if(xpi<xc+abs(raio(l)) && xpi>xc-abs(raio(l))) % check only if the x coordinate of the point is between 
                                                                                                                            % the values xc-radius and xc+radius
                                                                                                                            [teta_i,teta_f]=calcula_arco(x1,y1,x2,y2,xc,yc); % compute the arc between the line that defines the arc and
                                                                                                                            % the horizontal direction (-pi<teta<pi)
                                                                                                                            tetac(1)=acos((xpi-xc)/abs(raio(l))); % first intersection of the horizontal line that cross the 
                                                                                                                            % point with the circunference (angle between the radius that cross the intersection point and
                                                                                                                                                            % the horizontal direction)
                                                                                                                            tetac(2)=-tetac(1);  % angle of the second intersection 
                                                                                                                            y(1)=abs(raio(l))*sin(tetac(1))+yc; % y coordinate of the first intersection point
                                                                                                                            y(2)=abs(raio(l))*sin(tetac(2))+yc; % y coordinate of the second intersection point
                                                                                                                            for k=1:2 % check if the angles of the two intersection points are between the angle of the initial and the
                                                                                                                                % final angles of the arc. If yes so the vertical line from the the point to the botton direction
                                                                                                                                % intercept the arc (the counter should be increased).
                                                                                                                                if(raio(l)>0) % the center is on the left of the arc (from the initial point to the final point)
                                                                                                                                    if(teta_f>teta_i) 
                                                                                                                                        if(tetac(k)>teta_i && tetac(k)<teta_f)
                                                                                                                                            if(y(k)<ypi)
                                                                                                                                                interv=interv+1;
                                                                                                                                            end;
                                                                                                                                        end;
                                                                                                                                    else % teta_f<teta_i 
                                                                                                                                        if(tetac(k)>teta_i || tetac(k)<teta_f)
                                                                                                                                            if(y(k)<ypi)
                                                                                                                                                interv=interv+1;
                                                                                                                                            end;
                                                                                                                                        end
                                                                                                                                    end;
                                                                                                                                else % raio(l) < 0 the center is on the right of the arc (from the initial point to the final point)
                                                                                                                                    if(teta_i > teta_f)
                                                                                                                                        if(tetac(k)>teta_f && tetac(k)<teta_i)
                                                                                                                                            if(y(k)<ypi)
                                                                                                                                                interv=interv+1;
                                                                                                                                            end;
                                                                                                                                        end;
                                                                                                                                    else % teta_i < teta_f
                                                                                                                                        if(tetac(k)>teta_f || tetac(k)<teta_i)
                                                                                                                                            if(y(k)<ypi)
                                                                                                                                                interv=interv+1;
                                                                                                                                            end;
                                                                                                                                        end;
                                                                                                                                    end;
                                                                                                                                end;
end
end;
end;
l=l+1;
end;


if rem(interv,2) ~= 0	% Resto da divisão de interv por 2
    ponto = 'interno';
else
    ponto = 'externo';
end;function telles(gamm,eet)
  
    eest = eet^2 - 1;
    term1 = eet*eest + abs(eest);
    if term1 < 0
        term1 = (-term1)^(1/3);
        term1 = -term1;
    else
        term1 = term1^(1/3);
    end
  
    term2 = eet*eest - abs(eest);
    if term2 < 0
        term2 = (-term2)^(1/3);
        term2 = -term2;
    else
        term2 = term2^(1/3);
    end
    GAMM = term1 + term2 + eet;
  
  
    Q = 1 + 3*GAMM.^2;
    A = 1./Q;
    B = -3.*GAMM./Q;
    C = 3.*GAMM.^2./Q;
    D = -B;
  
    eta = A.*gamm.^3 + B.*gamm.^2 + C.*gamm + D;
    Jt = 3.*A.*gamm.^2 + 2.*B.*gamm + C;
    return eta,Jt
end
% Função para testar a proximidade de um ponto interno ao contorno
% Se o ponto estiver mais proximo que o desejado (d<d_min) o ponto
% é reprovado (aceita = não). Caso contrrário (aceita = sim)
%
%   Autor: Frederico Lourenço
%   Data de criação setembro de 1999 
%   Revisão 0.0

function [aceita] = aceita_dist(xpi,ypi,xl1,yl1,xl2,yl2,d_min,raio)
  
    pc = 0;					% Contador de pontos do contorno
        npc = length(xl1);	% N. de pontos do contorno
            aceita = 'sim';
          
            % Verificando a proximidade do ponto interno aos pontos do contorno		
                while (strcmp(aceita,'sim') && pc < npc)
                    pc = pc + 1;
                    x1 = xl1(pc);
                    y1 = yl1(pc);
                    d = sqrt((xpi-x1)^2+(ypi-y1)^2);
                    if d < d_min
                        aceita = 'não';
                    end;
                end;
              
                % Verificando a proximidade às linhas do contorno
                    l = 0;	% Contador das linhas do contorno
                      
                        while (strcmp(aceita,'sim') && l < npc)
                            l = l + 1;
                            x1 = xl1(l); y1 = yl1(l);
                            x2 = xl2(l); y2 = yl2(l);   
                            if(raio(l)==0) % The segment is a straight line
                              
                                if (xpi > x1 && xpi < x2) || (xpi > x2 && xpi < x1)
                                    m = (y2-y1)/(x2-x1);		% Angular coeficient of the recta
                                    yi = m*(xpi-x1)+y1;			% y coordinate of the intersection of a recta normal to the current boundary
                                    % recta that cross the point
                                    dy = ypi - yi; 
                                    d = abs(dy*cos(atan(m))); % distance from the point to the recta
                                    if d < d_min 
                                        aceita = 'não';    
                                    end;
                                end;
                              
                                if (ypi > y1 && ypi < y2) || (ypi > y2 && ypi < y1)
                                    if x1 == x2
                                        d = abs(xpi-x1);
                                    else 
                                        m = (y2-y1)/(x2-x1);			% Angular coeficient of the recta
                                        xi = 1/m*(ypi-y1)+x1;	% x coordinate of the intersection of a recta normal to the current boundary
                                        % recta that cross the point            
                                        dx = xpi - xi;
                                        d = abs(dx*sin(atan(m))); % distance from the point to the recta
                                    end;
                                    if d < d_min
                                        aceita = 'não';    
                                    end;
                                end;
                            else    % The segment is an arc
                                [xc,yc]=calcula_centro(x1,y1,x2,y2,raio(l)); % Center of the arc
                                [teta_i,teta_f]=calcula_arco(x1,y1,x2,y2,xc,yc); % angle of the lines that defines the arc with the horizontal direction
                                [teta_i,teta_p]=calcula_arco(x1,y1,xpi,ypi,xc,yc);  % teta_p angle of the line that cross the center point and the
                                % internal point with the horizontal direction
                                if(raio(l)>0) % The center is in the left side of the arc (from the initial to the end point)
                                    if(teta_f>teta_i) 
                                        if(teta_p>teta_i && teta_p<teta_f)
                                            d=abs(raio(l)-sqrt((xpi-xc)^2+(ypi-yc)^2)); % distance from the point to the arc
                                            if d < d_min
                                                aceita = 'não';    
                                            end;
                                        end;
                                    else % teta_f<teta_i
                                        if(teta_p>teta_i || teta_p<teta_f)
                                            d=abs(raio(l)-sqrt((xpi-xc)^2+(ypi-yc)^2)); % distance from the point to the arc
                                            if d < d_min
                                                aceita = 'não';    
                                            end;
                                        end;
                                    end;
                                else % raio(l) < 0 % The center is in the right side of the arc (from the initial to the end point)
                                    if(teta_i > teta_f)
                                        if(teta_p>teta_f && teta_p<teta_i)
                                            d=abs(abs(raio(l))-sqrt((xpi-xc)^2+(ypi-yc)^2)); % distance from the point to the arc
                                            if d < d_min
                                                aceita = 'não';    
                                            end;
                                        end;
                                    else % teta_i < teta_f
                                        if(teta_p>teta_f || teta_p<teta_i)
                                            d=abs(abs(raio(l))-sqrt((xpi-xc)^2+(ypi-yc)^2)); % distance from the point to the arc
                                            if d < d_min
                                                aceita = 'não';    
                                            end;
                                        end;
                                    end;
                                end;         
                            end;
                        end;


% Retorna a aprovação ou não do ponto interno
