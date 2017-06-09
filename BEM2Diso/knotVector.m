function [ knots ] = knotVector( ctrlPoints,grCurva )
%Calcula valores para o vetor nós
%   Matheus Augusto Correia - 04/08/15


    n=size(ctrlPoints,1);%número depontos de controle

    m=n+grCurva+1;%calcula o número de pontos
    
    %Calculo do vetor nós
    knots=zeros(1,m);
        
        for i=(m-grCurva):m
            knots(i)=1;
        end
       
        if rem(m,2)==0
            knots(m/2)=0.5;
            knots(m/2+1)=0.5;
            for i=(m/2+1):m-grCurva-2
                knots(i+1)=(knots(i)+1)/2;
            end
            for i=m/2:-1:grCurva+3
                knots(i-1)=knots(i)/2;
            end
        else
            knots(m/2+0.5)=0.5;
            for i=(m/2+0.5):m-grCurva-2
                knots(i+1)=(knots(i)+1)/2;
            end
            for i=(m/2+0.5):-1:grCurva+3
                knots(i-1)=knots(i)/2;
            end
        end
end

