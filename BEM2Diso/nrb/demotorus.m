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

function demotorus

% DEMOTORUS: A second demonstration of surface construction
% by revolution.
%

% D.M. Spink


sphere = nrbrevolve(nrbcirc(1,[],0.0,pi),[0.0 0.0 0.0],[1.0 0.0 0.0]);
nrbplot(sphere,[20 20],'light','on');
title('Ball and torus - surface construction by revolution');
hold on;
torus = nrbrevolve(nrbcirc(0.2,[0.9 1.0]),[0.0 0.0 0.0],[1.0 0.0 0.0]);
nrbplot(torus,[20 10],'light','on');
nrbplot(nrbtform(torus,vectrans([-1.8])),[20 10],'light','on');
hold off;

