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

function democurve

% Shows a simple test curve.
% 


crv = nrbtestcrv;

% plot the control points
plot(crv.coefs(1,:),crv.coefs(2,:),'ro');
title('Arbitrary Test 2D Curve.');
hold on;
plot(crv.coefs(1,:),crv.coefs(2,:),'r--');

% plot the nurbs curve
nrbplot(crv,48);
hold off;

crv.knots(4)=0.1;
figure
% plot the control points
plot(crv.coefs(1,:),crv.coefs(2,:),'ro');
title('Arbitrary Test 2D Curve.');
hold on;
plot(crv.coefs(1,:),crv.coefs(2,:),'r--');

% plot the nurbs curve
nrbplot(crv,48);
hold off;
