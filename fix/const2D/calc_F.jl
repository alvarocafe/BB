function calc_F(NOS,CW,FR,fc,finc)
# Evaluates the influence of any concentrated source or incident plane wave of the problem
nnos=size(NOS,1);	# Number of physical nodes
q = complex(zeros(nnos));	# Allocates the vector for the concentrated sources
inc = complex(zeros(nnos));	# Allocates the vector for incident plane waves
for i = 1:nnos	# Loop over the nodes
  xd=NOS[i,2]; # x coordinate of the source point
  yd=NOS[i,3]; # y coordinate of the source point
  q[i]=calc_q(xd,yd,fc,FR,CW);	# Evaluates the influence of the concentrated sources
  inc[i] = calc_inc(xd,yd,finc,FR,CW);	# Evaluates the influence of the incident plane wave
end
return q, inc
end
