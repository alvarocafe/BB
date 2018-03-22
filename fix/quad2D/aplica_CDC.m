function [A,b,T_PR]=aplica_CDC(G,H,CDC,ELEM)
% Aplica as condições de contorno (troca as colunas das matrizes G e H)

n_el=size(CDC,1);
todos_valores=zeros(1,3*n_el);
T_PR=zeros(1,5);
cont=0;

% Cria a variável T_PR que contém os nós onde a temperatura é prescrita
% (conhecida).
% T_PR tem 5 colunas e o número de linhas é igual ao número de nós onde a
% temperatura é conhecida.
% T_PR=[a1,a2,a3,a4,a5]
% a1=número do nó com temperatura prescrita.
% a2=número do primeiro elemento com temperatura prescrita ao qual este nó
% pertence.
% a3=número local do nó neste elemento.
% a4: caso a temperatura também seja prescrita no segundo elemento a que
% este nó pertence, então, a4 conterá o número deste elemento, caso
% contrário, conterá zero.
% a5: caso a4 seja diferente de zero, a5 conterá o número local do nó no
% segundo elemento, caso contrário, conterá zero.

for el=1:n_el % for sobre os elementos para criar T_PR
    for no=1:3 % for sobre os nós do elemento el
        no_global=ELEM(el,no+1); % número global do nó
        tipoCDC=CDC(el,2*no);  % tipo da condição de contorno
        valorCDC=CDC(el,2*no+1); % valor da condição de contorno
        todos_valores(3*el-3+no)=valorCDC; % armazerna o valor da condição
                   % de contorno no vetor todos_valores
        if(tipoCDC==0) % se a temperatura é conhecida
            compartilha=0; % por enquanto não se sabe se a tempertura é 
            % também conhecida no segundo elemento a que este nó pertence
            i=1;
            while (~compartilha&&i<=cont&&cont>0) % verifica se o nó global
                % já está presente em T_PR. Quando ele encontra,
                % compartilha se torna igual a 1 e o while pára.
                if(no_global==T_PR(i,1)) % Se sim, o nó global já está 
                    % presente em T_PR. As colunas 4 e 5 são preenchidas e
                    % compartilha se torna igual a 1.
                    compartilha=1;
                    T_PR(i,4)=el;
                    T_PR(i,5)=no;
                end
                i=i+1;
            end
            if(~compartilha) % Se compartilha continua zero, então o nó 
                % ainda não foi inserido em T_PR. Neste caso, as três
                % primeiras colunas de T_PR são preenchidas.
                cont=cont+1;
                T_PR(cont,1)=no_global;
                T_PR(cont,2)=el;
                T_PR(cont,3)=no;
            end
        end
    end
end

n_temp_pr = length(T_PR(:,1)); % Número de nós com temperatura conhecidas
for i=1 : n_temp_pr % for sobre os nós com temperatura conhecida
    i_no=T_PR(i,1); % número do nó com temperatura conhecida
    i_el=T_PR(i,2); % primeiro elemento que contém o nó
    i_no_loc=T_PR(i,3); % número local do nó neste elemento
    ind_H=i_no; % índice da coluna da matriz H que será trocada
    ind_G=3*i_el+i_no_loc-3; % índice da coluna da matriz G que será 
                 %  trocada
    troca = G(:,ind_G); % armazena a coluna de G que será trocada
    G(:,ind_G) = -H(:,ind_H); % substitui a coluna de H na de G
    H(:,ind_H) = -troca; % Substitui a de G na de H
    if(T_PR(i,4)~=0) % Se a temperatura também é conhecida no segundo 
         %   elemento
        i_el=T_PR(i,4); % Número do segundo elemento 
        i_no_loc=T_PR(i,5); % número local deste nó no segundo elemento
        ind_G=3*i_el+i_no_loc-3; % índice da coluna G que será atribuído 
                  % zero
        H(:,ind_H)=H(:,ind_H)-G(:,ind_G); % soma o valor da coluna G na 
                  % coluna H
        G(:,ind_G)=0; % atribui zero na coluna G
    end;
end;
b=G*todos_valores'; % Cálculo do vetor b
A=H;
return
