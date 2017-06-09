function [H,G]=calcula_PropGeom(n,collocCoord,collocPts,nnos,crv,dcrv, ...
    npgauss,kmat)

H=zeros(length(collocCoord));
G=zeros(length(collocCoord));
% P=zeros(length(collocCoord));
[qsi1,w1]=Gauss_Legendre(-1,1,npgauss); % Calcula pesos e pontos de Gauss
% usados na integra��o de I= int F n.r/r dGama


for i=1:n % Loop sobre o n�mero de curvas
    n_el = size(crv(i).elRange,1);	% N�mero total de elementos
    % qsi1 e qsi2 = pontos de Gauss
    % w1 e w2 = pesos de Gauss
    
    
    for j = 1 : n_el	% Percorre os elementos
        for k=1:length(collocCoord) % Loop sobre os pontos de coloca��o
            xfonte=collocCoord(k,1);
            yfonte=collocCoord(k,2);
            % Integrais ao longo dos elementos (que agora se sobrep�e)
            [g,h,id]=integra_elem(xfonte,yfonte,crv(i),dcrv(i),crv(i).elRange(j,:),...
                qsi1,w1,kmat,collocPts(k)); % Integra��o sobre o
            %  element
            H(k,id+1+nnos(i))=H(k,id+1+nnos(i))+h;
            G(k,id+1+nnos(i))=G(k,id+1+nnos(i))+g;
        end
    end
end
return
