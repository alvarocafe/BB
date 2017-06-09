function contorno=calc_ncont(SEGMENTOS)
num_seg=size(SEGMENTOS,1);
t=1;
p2=0;
p_ini = SEGMENTOS(1,2);
icont=1; % Índice de contornos
contorno(1,1)=p_ini;
while(t<num_seg)  	% While over all lines
    while(p2~=p_ini)
        p1  = SEGMENTOS(t,2);
        p2  = SEGMENTOS(t,3);
        if p2 == p_ini
            if t < num_seg
                p_ini = SEGMENTOS(t+1,2);
                icont=icont+1;
                contorno(icont,1)=p_ini;
                contorno(icont-1,2)=p_ini-contorno(icont-1,1);
            end;
        end;
        t=t+1;
    end;                                  %end of while p2
end;
if(icont>1)
    contorno(icont,2)=num_seg-sum(contorno(1:icont-1,2));
else
    contorno(1,2)=num_seg;
end
return