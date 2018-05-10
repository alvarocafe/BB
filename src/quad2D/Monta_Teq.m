function [T,q,q_vet]=Monta_Teq(CDC,ELEM,x)

nelem=length(ELEM(:,1)); %Número de elementos
tipoCDCultimo=CDC(nelem,2);
n_pontos=length(x);
for el = 1:nelem  % Corre os elementos
    no1 = ELEM(el,2); % Número do primeiro nó do elemento
    no2 = ELEM(el,3); % Número do segundo nó do elemento
    no3 = ELEM(el,4); % Número do segundo nó do elemento
    tipoCDC=CDC(el,2); % Tipo da condição de contorno no elemento ij
    valorCDC=CDC(el,3); % Valor da condição de contorno no elemento ij
    if tipoCDC==0   % Temperatura conhecida
        T([no1,no2,no3])=[valorCDC valorCDC valorCDC]; % Atribui os valores
        % das condições de contorno ao vetor T
        q(el,1:4)=[el x(no1) x(no2) x(no3)]; % Atribui os valores calculados
        q_vet(3*el-2:3*el)=[x(no1) x(no2) x(no3)];
        % ao vetor q
    else % fluxo conhecido
        q(el,1:4)=[el valorCDC valorCDC valorCDC]; % Atribui os valores
        % das condições de contorno ao vetor q
        q_vet(3*el-2:3*el)=[valorCDC valorCDC valorCDC];
        if el==1 && tipoCDCultimo==1
            T([no1 no2])=x([no1 no2 no3]);
        elseif(el~=nelem)
            T([no2 no3])=[x(no2) x(no3)];
        else
            T(no2)=x(no2);
        end
    end
end
T(2*nelem+1:n_pontos)=x(2*nelem+1:n_pontos);