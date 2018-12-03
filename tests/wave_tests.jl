## Test cases for wave problems
# This script tests the BEM models against known problems
#for which analytical solutions are readily obtained.
#################### Test case n ####################
### Problem description
### Analytical solution
### BEM model
### return error
include("../src/const2D/const2D.jl")
include("../src/nurbs2D/nurbs2D.jl")
using .const2D, .nurbs2D, SpecialFunctions

println("Running tests for wave problems...")

#################### Test case 1 ####################
### Acoustic square
# Consider a square acoustic domain in which the speed of
#sound is approximately 344 [m/s] inside the domain and
#there are two opposite walls which are rigid and two
#open-ended, such that one of them is excited at a frequency
# ω [rads], so that the wave number is k = ω/c, with
#relative amplitude P = 1.
### Analytical solution
phi_square(x) = sin.(k.*x)
q_square(x) = -k.*cos(k.*x)
### BEM model
function square(ne=10,L=1,k=1)
    #L = 1; # length of the square
    #k = 1; # wave number of the problem
    #The points and segments which describe this geometry are
    POINTS = [1 0 0
	      2 L 0
	      3 L L
	      4 0 L];
    SEGMENTS = [1 1 2 0
	        2 2 3 0
	        3 3 4 0
	        4 4 1 0];
    # Each segment will be meshed by ne elements
    #ne = 100;
    MESH = [1 ne
	    2 ne
	    3 ne
	    4 ne];
    BCSeg = [1 1 0
	     2 0 0
	     3 1 0
	     4 0 1];
    PONTOS_dom = [1 L/2 L/2]
    fc = [0 0 0]
    # Conventional method with full influence matrices
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    info = [NOS_GEO,NOS,ELEM,CDC]
    tic()
    phi, q, phi_dom, phi_dom = const2D.solve(info,PONTOS_dom,fc,BCSeg,k)
    t = toq()
    # Conventional method with approximated influence matrices
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    info = [NOS_GEO,NOS,ELEM,CDC]
    tic()
    phiH,qH,phi_domH,phi_domH = const2D.solveH(info,PONTOS_dom,fc,BCSeg,k)
    tH = toq()
    return t, tH, phi, phiH, phi_dom
end

#################### Test case 2 ####################
### Vibrating cylinder
# An acoustically rigid cylinder is let to vibrate in an
#infinite acoustic domain at frequency ω so that the wavenumber is k = ω/c.
### Analytical solution
#The analytical solution for the vibratin cylinder obtained
#by applying variable separation in cylindrical coordinates
#of the Helmholtz equation is
phi_cylinder(k,r,x) = (1/k).*(besselh(0,2,k.*x)./besselh(1,2,k.*r));
#at a distance x from the cylinder of radius r.
### BEM model
function cylinder(ne = 100,r=0.5,c=[0 0],k=1)
    t=0; tH=0; tiso=0; ϵ=0; ϵH=0; ϵiso=0; phi=0; q=0; phi_dom=0;  phiH=0; qH=0; phi_domH=0;  phiiso=0; qiso=0; domiso=0;
    POINTS =[1  c[1,1]-r	c[1,2]
       	     2	c[1,1]		c[1,2]+r
       	     3	c[1,1]+r	c[1,2]
       	     4	c[1,1]		c[1,2]-r];
    SEGMENTS = [1 1 2 -r
                2 2 3 -r
                3 3 4 -r
                4 4 1 -r];
    BCSeg = [1 1 1 0
             2 1 1 0
             3 1 1 0
             4 1 1 0];
    MESH = [1 ne
            2 ne
            3 ne
            4 ne];
    x = 10*r
    PONTOS_dom = [1 0 x]
    fc = [0 0 0]
    # # Conventional method with full influence matrices
     NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
     nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
     info = [NOS_GEO,NOS,ELEM,CDC]
    # tic()
    # phi, q, phi_dom, phi_dom = const2D.solve(info,PONTOS_dom,fc,BCSeg,k)
    # t = toq()
    # ϵ = abs.(sqrt.(((phi_dom .- phi_cylinder(k,r,x)).^2)./phi_cylinder(k,r,x).^2))
    # Conventional method with approximated influence matrices
    tic()
    phiH,qH,phi_domH,phi_domH = const2D.solveH(info,PONTOS_dom,fc,BCSeg,k)
    tH = toq()
    ϵH = abs.(sqrt.(((phi_domH .- phi_cylinder(1,0.5,5)).^2)./phi_cylinder(1,0.5,5).^2))
    # # Isogeometric BEM with full influence matrices
    # crv = nurbs2D.format_dad(POINTS,SEGMENTS,MESH)
    # dcrv=map(x->nurbs2D.nrbderiv(x),crv)
    # n = length(crv)
    # z=0;
    # for k=1:n
    #     for i=1:crv[k].number
    #         z=z+1
    #     end
    # end
    # numcurva=zeros(Integer,z)
    # collocPts=zeros(z)
    # CDC=zeros(z,3)
    # collocCoord=zeros(z,3)
    # z=0;
    # nnos=zeros(Integer,n)
    # for k=1:n
    #     p=crv[k].order-1;
    #     nnos[k]=crv[k].number;
    #     valorCDC=BCSeg[k,3];
    #     tipoCDC=BCSeg[k,2];
    #     for i=1:crv[k].number
    #         z=z+1;
    #         numcurva[z]=k;
    #         collocPts[z]=sum(crv[k].knots[(i+1):(i+p)])/p;
    #         if(i==2)
    #             collocPts[z-1]=(collocPts[z]+collocPts[z-1])/2;
    #         end
    #         if(i==nnos[k])
    #             collocPts[z]=(collocPts[z]+collocPts[z-1])/2;
    #         end

    #         CDC[z,:] = [z,tipoCDC,valorCDC];
    #     end
    # end
    # nnos2=cumsum([0 nnos'],2)
    # E=zeros(length(collocPts),length(collocPts));
    # for i=1:length(collocPts)
    #     collocCoord[i,:]=nurbs2D.nrbeval(crv[numcurva[i]], collocPts[i]);
    #     B, id = nurbs2D.nrbbasisfun(crv[numcurva[i]],collocPts[i])
    #     E[i,id+nnos2[numcurva[i]]]=B
    # end
    # info = collocCoord,nnos2,crv,dcrv,E
    # tic()
    # phiiso, qiso, domiso, domiso = nurbs2D.solve(info,PONTOS_dom,fc,CDC,k)
    # tiso = toq()
    # ϵiso = abs.(sqrt.(((domiso .- phi_cylinder(k,r,x)).^2)./phi_cylinder(k,r,x).^2))
    return t, tH, tiso, ϵ, ϵH, ϵiso, phi, q, phi_dom,  phiH, qH, phi_domH,  phiiso, qiso, domiso
end
