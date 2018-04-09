function [dTdx,dTdy]=calc_dT_pint(ELEM,NOS,PONTOS_INT,T,q,k,fc)
% Calcula a temperatura nos pontos internos
n_pint=length(PONTOS_INT(:,1)); % Numero de pontos internos
n_elem=length(T); % Numero de elementos
Sx_int=zeros(n_pint,n_elem);
Sy_int=zeros(n_pint,n_elem);
Dx_int=zeros(n_pint,2*n_elem);
Dy_int=zeros(n_pint,2*n_elem);
dqdx_fc=zeros(n_pint,1);
dqdy_fc=zeros(n_pint,1);

npg=16;
[qsi,w]=Gauss_Legendre(-1,1,npg);


for i=1:n_pint % Laço sobre os pontos internos
    x_fonte=PONTOS_INT(i,2); % Coordenada x do ponto fonte
    y_fonte=PONTOS_INT(i,3); % Coordenada y do ponto fonte
    for j=1:n_elem
        x1=NOS(ELEM(j,2),2); % Coordenada x do inicio do elemento
        y1=NOS(ELEM(j,2),3); % Coordenada y do inicio do elemento
        x2=NOS(ELEM(j,3),2); % Coordenada x do final do elemento
        y2=NOS(ELEM(j,3),3); % Coordenada y do final do elemento
        
        [dx,dy,sx,sy]=calc_DeS(x1,y1,x2,y2,x_fonte,y_fonte,k,qsi,w); % Chama a
        % functio para cálculo de H e G quando o
        % ponto fonte nao pertence ao elemento
        for nolocal = 1 : 2
            noglobal = ELEM(j,nolocal+1); %Índice da matriz global H
            Sx_int(i,noglobal) = Sx_int(i,noglobal) + sx(nolocal);
            Sy_int(i,noglobal) = Sy_int(i,noglobal) + sy(nolocal);
        end;
        Dx_int(i,2*j-1:2*j) = dx;
        Dy_int(i,2*j-1:2*j) = dy;
    end
    [dqdx_fc(i),dqdy_fc(i)]=calc_dq(x_fonte,y_fonte,fc,k);    
end

dTdx=Dx_int*q'-Sx_int*T-dqdx_fc; % Vetor que contem a temperatura nos pontos internos
dTdy=Dy_int*q'-Sy_int*T-dqdy_fc; % Vetor que contem a temperatura nos pontos internos