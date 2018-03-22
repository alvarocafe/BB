function calc_inc(xd,yd,finc,k)
# Evaluates the influence of the incident plane waves described by finc
# finc = [flag A x y], where flag is either 1 or 0 and indicates the presence of incident plane waves, A is the amplitude of the wave, x and y are the direction of propagation of the wave (sqrt(x^2 + y^2) = 1)
n_inc = size(finc,1);	# Number of incident waves
inc = complex(0,0);	# Inlfuence of all incident plane waves
for i = 1:n_inc	# Loop over the incident waves
  inc += finc[i,4]*exp(complex(0,k*(finc[i,2]*xd + finc[i,3]*yd)));	# Sums the influence of each plane wave
end
return inc
end
