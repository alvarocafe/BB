function  [g,h,id]=integra_elem(xf,yf,crv,dcrv,range,...
    qsi1,w1,k,cpt)

% Integra��o sobre os elementos (integral I)

% N�mero de pontos de Gauss usados na integra��o de I
npg = length(qsi1);

g = zeros(1,crv.order); % Inicializa o per�metro
h = zeros(1,crv.order); % Inicializa o per�metro
dudqsi=(range(2)-range(1)) / 2;
eet = (range(1)+range(2)-2*cpt)/(range(1)-range(2));
[eta,Jt]=telles(qsi1,eet);
for i = 1 : npg % Percorre os pontos de integra��o
    qsi_param=convertToParamSpace(eta(i), range);
    [B, id] = nrbbasisfun (qsi_param,crv);    
    [pt,dxydu] = nrbdeval(crv, dcrv, qsi_param);
    x=pt(1,:);y=pt(2,:);
    dgamadu=norm(dxydu);    
    sx=dxydu(1)/dgamadu; % Component x do vetor tangente
    sy=dxydu(2)/dgamadu; % Componente y do vetor tangente
    nx=sy; % Componente x do vetor normal unit�rio
    ny=-sx; % Componente y do vetor normal unit�rio
    [Tast,qast]=calc_solfund(xf,yf,x,y,nx,ny,k);
    % C�lculo das integrais
    h = h + B*qast* dgamadu *dudqsi* w1(i)* Jt(i); % Integral da
    g = g + B*Tast* dgamadu *dudqsi* w1(i)* Jt(i);  % Integral da
end


%Valor de qsi no ponto fonte (pode ser fora do elemento se o ponto fonte
% n�o pertencer ao elemento. Ver o artigo da Transformada de Telles.
% if abs(x1-x2)~=0 % >= 1e-6
%     eet = (x1+x2-2*x_f)/(x1-x2);
% else
%     eet = (y1+y2-2*y_f)/(y1-y2);
% end
% for ii = 1 : npg
%     [PT,JT] = telles(PG(ii),eet);
%     x_campo = x_med + PT*d_x/2;
%     y_campo = y_med + PT*d_y/2;
%     [FHe,FGe,FCij] = sol_fund_pla(x_f,y_f,x_campo,y_campo,nx,ny,Material);
%     He = FHe*JT*WG(ii)*Jac + He;
%     Cije = FCij*JT*WG(ii)*Jac + Cije;
%     Ge = FGe*JT*WG(ii)*Jac + Ge;
% end

return
