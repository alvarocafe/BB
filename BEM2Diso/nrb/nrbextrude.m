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

function srf = nrbextrude(curve,vector)

%
% NRBEXTRUDE: Construct a NURBS surface by extruding a NURBS curve.
% 
% Calling Sequence:
% 
%   srf = nrbextrude(crv,vec);
% 
% Parameters:
% 
%   crv		: NURBS curve to extrude, see nrbmak.
% 
%   vec		: Vector along which the curve is extruded.
% 
%   srf		: NURBS surface constructed.
% 
% Description:
% 
%   Constructs a NURBS surface by extruding a NURBS curve along a defined 
%   vector. The NURBS curve forms the U direction of the surface edge, and
%   extruded along the vector in the V direction. Note NURBS surfaces cannot
%   be extruded.
% 
% Examples:
% 
%   Form a hollow cylinder by extruding a circle along the z-axis.
%
%   srf = nrbextrude(nrbcirc, [0,0,1]);

%  D.M. Spink
%  Copyright (c) 2000.

if iscell(curve.knots)
  error('Nurb surfaces cannot be extruded!');
end

if nargin < 2
  error('Error too few input arguments!');
end

if nargin == 3
  dz = 0.0;
end

coefs = cat(3,curve.coefs,vectrans(vector)*curve.coefs);
srf = nrbmak(coefs,{curve.knots, [0 0 1 1]});

