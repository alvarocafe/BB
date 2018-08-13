## Test cases for heat conducticity problems
# This script tests the BEM models against
#known problems for which analytical solutions
#are readily obtained.
########## Test case number ##########
### Problem description
### Analytical solution
### BEM model
### return error


########## Test case 1 ##########
# Let a square domain of length L [m] and heat
#conductivity k [W/(m.K)], initially at
#uniform constant temperature be subjected to
#a temperature gradient on its boundaries on
#opposites sides in the x direction.
L = 1; # length of the square
k = 1; # heat conductivity of the material
#The points and segments which describe this
#geometry are
POINTS = [1 0 0
	  2 L 0
	  3 L L
	  4 0 L];
SEGMENTS = [1 1 2 0
	    2 2 3 0
	    3 3 4 0
	    4 4 1 0];
# Each segment will be meshed by ne elements
ne = 1;
MESH = [1 ne
	2 ne
	3 ne
	4 ne];
# BCSeg =[NE BT V] are the boundary conditions
#at each segment. [NE BT V] are
#NE: element number
#BT: boundary condition type
# BT = 0 if the temperature is known
# BT = 1 if the temperature gradient is known
#V: value of boundary variable
BCSeg = [1 1 0
	 2 0 0
	 3 1 0
	 4 0 1];
NPX = 11; # Number of domain points on x and y
NPY = 11;
### Analytical solution
### BEM modules test
# Each module will be tested in this function.

function test_modules(ϵ=10^(-6))
    ## TODO: automate module search and testing.
    # Constant elements
    include("../src/const2D/const2D.jl")
    using const2D
    NOS_GEO,NOS,ELEM,CDC,normal = const2D.format_dad(POINTS,SEGMENTS,MESH,BCSeg) # Apply the discretization technique and builds the problems matrices for the geometrical points, physical nodes, elements' connectivity and boundary conditions
    nnos = size(NOS,1)  # Number of physical nodes, same as elements when using constant elements
    G,H=const2D.cal_GeHpot(NOS,NOS_GEO,ELEM,k,fc,qsi,w);
    A,b = const2D.aplica_CDC(G,H,CDC);
    x = A\b # Solves the linear system
    T,q = const2D.monta_phieq(CDC,x) # Applies the boundary conditions to return the velocity potential and flux

    
    return const2D.error
end
