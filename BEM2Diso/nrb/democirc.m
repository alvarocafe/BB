% function democirc
% Demonstration of a circle and arcs in the x-y plane.
%

% D.M. Spink
% Copyright (c) 2000

crv = nrbcirc(1,[],deg2rad(45),deg2rad(315));
dcrv=nrbderiv(crv);
ddcrv=nrbderiv(dcrv);
us=0:0.1:1;
ps=nrbeval(crv,us);
P=[-sqrt(2)/2; sqrt(2)/2;0];
[~,ind]=min(sqrt(sum((ps-P).^2,1)));
u=us(ind);
for i=1:1000
    corr=-(nrbdeval(crv,dcrv,u)'*(nrbeval(crv,u)-P))/(nrbdeval(dcrv,ddcrv,u)'*(nrbeval(crv,u)-P)+norm(nrbdeval(crv,dcrv,u))^2);
    u=u+corr;
    if abs(corr)<1e-6
        break
    end
end