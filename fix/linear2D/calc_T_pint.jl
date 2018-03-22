function T_pint=calc_T_pint(ELEM,NOS,PONTOS_INT,T,q,k,fc)
% Calcula a temperatura nos pontos internos
n_pint=length(PONTOS_INT(:,1)); % Numero de pontos internos
n_elem=length(T); % Numero de elementos
H_int=zeros(n_pint,n_elem);
G_int=zeros(n_pint,2*n_elem);
q_fc=zeros(n_pint,1);

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
        [g,h]=calc_gh_nsing(x1,y1,x2,y2,x_fonte,y_fonte,k,qsi,w); % Chama a
        % functio para cálculo de H e G quando o
        % ponto fonte nao pertence ao elemento
        for nolocal = 1 : 2
            noglobal = ELEM(j,nolocal+1); %Índice da matriz global H
            H_int(i,noglobal) = H_int(i,noglobal) + h(nolocal);
        end;
        G_int(i,2*j-1:2*j) = g;
    end
    q_fc(i)=calc_q(x_fonte,y_fonte,fc,k);
end

T_pint=H_int*T-G_int*q'+q_fc; % Vetor que contem a temperatura nos pontos internos