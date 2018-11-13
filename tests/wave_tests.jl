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
using const2D

println("Running tests for BEM_base...")

#################### Test case 1 ####################
### Problem description
# Consider a square acoustic domain. The speed of
#sound is approximately 344 [m/s] inside the domain and
#there are two opposite walls which are rigid and two
#open-ended, such that one of them is excited at a frequency
# ω [rads], so that the wave number is k = ω/c, with
#relative amplitude P = 1.
### Analytical solution
phi_an(x) = sin.(k.*x)
q_an(x) = -k.*cos(k.*x)
### BEM model
L = 10; # length of the square
k = 1; # wave number of the problem
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
ne = 100;
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



function test_wave()
    nei = 50; nef = 300; passo = 10; nev = zeros(round(Int,(nef-nei)/passo + 1));
    t = zeros(round(Int,(nef-nei)/passo +1))
    i=0
    for ne = nei:passo:nef
        i+=1
        nev[i] = ne
        tic()
        MESH = [1 ne
	        2 ne
	        3 ne
	        4 ne];
        NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
        nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
        info = [NOS_GEO,NOS,ELEM,CDC]
        phi, q, phi_dom, phi_dom = const2D.solve(info,PONTOS_dom,fc,BCSeg,k)
        t[i] = toq()
    end
    tH = zeros(round(Int,(nef-nei)/passo +1))
    i=0
    for ne = nei:passo:nef
        i+=1
        nev[i] = ne
        tic()
        MESH = [1 ne
	        2 ne
	        3 ne
	        4 ne];
        NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg)
        nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
        info = [NOS_GEO,NOS,ELEM,CDC]
        tic()
        phiH,qH,phi_domH,phi_domH = const2D.solveH(info,PONTOS_dom,fc,BCSeg,k)
        tH[i] = toq()
    end
    return t, tH, error, errorH
end
# Run test function
#t, tH, error, errorH = test_wave();
