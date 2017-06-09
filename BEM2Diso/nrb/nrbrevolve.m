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

function surf = nrbrevolve(curve,pnt,vec,theta)

% 
% NRBREVOLVE: Construct a NURBS surface by revolving a NURBS curve.
% 
% Calling Sequence:
% 
%   srf = nrbrevolve(crv,pnt,vec[,ang])
% 
% Parameters:
% 
%   crv		: NURBS curve to revolve, see nrbmak.
% 
%   pnt		: Coordinate of the point used to define the axis
%               of rotation.
% 
%   vec		: Vector defining the direction of the rotation axis.
% 
%   ang		: Angle to revolve the curve, default 2*pi
% 
% Description:
% 
%   Construct a NURBS surface by revolving the profile NURBS curve around
%   an axis defined by a point and vector.
% 
% Examples:
% 
%   Construct a sphere by rotating a semicircle around a x-axis.
%
%   crv = nrbcirc(1.0,[0 0 0],0,pi);
%   srf = nrbrevolve(crv,[0 0 0],[1 0 0]);
%   nrbplot(srf,[20 20]);
%
% NOTE:
%
%   The algorithm:
%
%     1) vectrans the point to the origin (0,0,0)
%     2) rotate the vector into alignment with the z-axis
%
%     for each control point along the curve
%
%     3) determine the radius and angle of control
%        point to the z-axis
%     4) construct a circular arc in the x-y plane with 
%        this radius and start angle and sweep angle theta 
%     5) combine the arc and profile, coefs and weights.
%  
%     next control point
%
%     6) rotate and vectrans the surface back into position
%        by reversing 1 and 2.
%

%  D.M. Spink
%  Copyright (c) 2000.

if nargin < 3
  error('Not enough arguments to construct revolved surface');
end

if nargin < 4
  theta = 2.0*pi;
end

% Translate curve the center point to the origin
if isempty(pnt)
  pnt = zeros(3,1);
end

if length(pnt) ~= 3
  error('All point and vector coordinates must be 3D');
end

% Translate and rotate the curve into alignment with the z-axis
T  = vectrans(-pnt);
angx = vecangle(vec(1),vec(3));
RY = vecroty(-angx);
vectmp = RY*[vecnorm(vec(:));1.0];
angy = vecangle(vectmp(2),vectmp(3));
RX = vecrotx(angy);
curve = nrbtform(curve,RX*RY*T);

% Construct an arc 
arc = nrbcirc(1.0,[],0.0,theta);

% Construct the surface
coefs = zeros(4,arc.number,curve.number);
angle = vecangle(curve.coefs(2,:),curve.coefs(1,:));
radius = vecmag(curve.coefs(1:2,:));
for i = 1:curve.number  
  coefs(:,:,i) = vecrotz(angle(i))*vectrans([0.0 0.0 curve.coefs(3,i)])*...
          vecscale([radius(i) radius(i)])*arc.coefs;
  coefs(4,:,i) = coefs(4,:,i)*curve.coefs(4,i);
end
surf = nrbmak(coefs,{arc.knots, curve.knots});

% Rotate and vectrans the surface back into position
T = vectrans(pnt);
RX = vecrotx(-angy);
RY = vecroty(angx);
surf = nrbtform(surf,T*RY*RX);  




