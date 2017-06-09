%% Copyright (C) 2003 Mark Spink
%% 
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%% 
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%% 
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, see <http://www.gnu.org/licenses/>.

function demodercrv

% Demonstrates the construction of a general
% curve and determine of the derivative.
%




% make and draw nurbs test curve
crv = nrbtestcrv;
nrbplot(crv,48);
title('First derivatives along a test curve.');

npts = 9;
tt = linspace(0.0,1.0,npts);

dcrv = nrbderiv(crv);

% first derivative
[p1, dp] = nrbdeval(crv,dcrv,tt);

p2 = vecnorm(dp);

hold on;
plot(p1(1,:),p1(2,:),'ro');
h = quiver(p1(1,:),p1(2,:),p2(1,:),p2(2,:),0);
set(h,'Color','black');
hold off;




