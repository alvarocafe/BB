## Test cases for wave problems
# This script tests the BEM models against known problems
#for which analytical solutions are readily obtained.
#################### Test case n ####################
### Problem description
### Analytical solution
### BEM model
### return error
include("../src/const2D/const2D.jl")
# include("../src/nurbs2D/nurbs2D.jl")
# using nurbs2D
using .const2D, SpecialFunctions

println("Running tests for BEM_base...")

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
    return t, tH, phi, phiH
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
function cylinder(ne = 1000,r=0.5,c=[0 0],k=1)
    #r = 0.5; # cylinder radius
    #c = [0 0]; # cylinder center
    #ne = 100;
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
    PONTOS_dom = [1 0 10*r]
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
    return t, tH, phi_dom, phi_domH
end
