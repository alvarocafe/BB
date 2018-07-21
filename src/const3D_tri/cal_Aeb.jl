function cal_Aeb(b1,b2,arg)
    NOS, NOS_GEO, ELEM, k, CDC,qsi,w,qsi_tri,w_tri = arg;
    nnos = length(b1);
    nelem = length(b2)
    G=zeros(nnos,nelem); 
    H=zeros(nnos,nelem);     
    A = complex(zeros(nnos,nelem));
    b = complex(zeros(nnos, 1));
    ci=0
    for i in b1
        ci+=1
        xd=NOS[i,2]; 
        yd=NOS[i,3]; 
        zd=NOS[i,4];
        
        cj=0
        for j in b2
            cj+=1
	    tipoCDC = CDC[j,2]; 
	    
	    valorCDC = CDC[j,3];
	    nos = ELEM[j,2:4];		

            no1=ELEM[j,2]; 
            no2=ELEM[j,3]; 
            no3=ELEM[j,4]; 

            x1=NOS_GEO[no1,2]; 
            y1=NOS_GEO[no1,3]; 
            z1=NOS_GEO[no1,4]; 

            x2=NOS_GEO[no2,2]; 
            y2=NOS_GEO[no2,3]; 
            z2=NOS_GEO[no2,4]; 

            x3=NOS_GEO[no3,2]; 
            y3=NOS_GEO[no3,3]; 
            z3=NOS_GEO[no3,4]; 

            n = calc_vetnormal(x1,y1,z1,x2,y2,z2,x3,y3,z3); 
            
            if i==j 
                g,h=calcula_HeGs_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,qsi,w,k); 
            else 
                g,h=calcula_HeGns_POT(x1,y1,z1,x2,y2,z2,x3,y3,z3,xd,yd,zd,n,qsi_tri,w_tri,k); 
            end            
            if tipoCDC == 0 
        	A[ci,cj] = -g; 	
                if valorCDC == 0
                    b[ci,1] = b[ci,1] - 0; 
                else
                    b[ci,1] = b[ci,1] - h*valorCDC; 
                end                
            else
                A[ci,cj] = +h; 	
                if valorCDC == 0
                    b[ci,1] = b[ci,1] + 0; 
                else
                    b[ci,1] = b[ci,1] + g*valorCDC; 
                end
            end
        end
    end
    return A,b
end
