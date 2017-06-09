%% Copyright (C) 2009 Carlo de Falco
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

function idx = nrbnumbasisfun (points, nrb)

% NRBNUMBASISFUN:  Numbering of basis functions for NURBS
%
% Calling Sequence:
% 
%   N      = nrbbasisfun (u, crv)
%   N      = nrbbasisfun ({u, v}, srf)
%   N      = nrbbasisfun (p, srf)
%
%    INPUT:
%   
%      u or p(1,:,:)  - parametric points along u direction
%      v or p(2,:,:)  - parametric points along v direction
%      crv - NURBS curve
%      srf - NURBS surface
%   
%    OUTPUT:
%
%      N - Indices of the basis functions that are nonvanishing at each
%          point. size(N) == size(B)
%   

  if (   (nargin<2) ...
      || (nargout>1) ...
      || (~isstruct(nrb)) ...
      || (iscell(points) && ~iscell(nrb.knots)) ...
      || (~iscell(points) && iscell(nrb.knots) && (size(points,1)~=2)) ...
      )
    print_usage();
  end


  if (~iscell(nrb.knots))    %% NURBS curve
    
    iv  = findspan (nrb.number-1, nrb.order-1, points, nrb.knots);
    idx = numbasisfun (iv, points, nrb.order-1, nrb.knots);
    
  else                       %% NURBS surface

    if (iscell(points))
      [v, u] = meshgrid(points{2}, points{1});
      p(1,:,:) = u;
      p(2,:,:) = v;
      p = reshape(p, 2, []);
    else
      p = points;
    end
    
    idx = nrb_srf_numbasisfun__ (p, nrb); 

  end
  
end

function idx = nrb_srf_numbasisfun__ (points, nrb)
  
  m   = nrb.number(1)-1;
  n   = nrb.number(2)-1;
  
  npt = size(points,2);
  u   = points(1,:);
  v   = points(2,:);

  U   = nrb.knots{1};
  V   = nrb.knots{2};

  p   = nrb.order(1)-1;
  q   = nrb.order(2)-1;

  spu = findspan (m, p, u, U); 
  Ik  = numbasisfun (spu, u, p, U);

  spv = findspan (n, q, v, V);
  Jk  = numbasisfun (spv, v, q, V);
  
  for k=1:npt
    [Jkb, Ika] = meshgrid(Jk(k, :), Ik(k, :)); 
    idx(k, :)  = sub2ind([m+1, n+1], Ika(:)+1, Jkb(:)+1);
  end

end

%!test
%! p = 2;
%! q = 3;
%! mcp = 2; ncp = 3;
%! knots = {[zeros(1,p), linspace(0,1,mcp-p+2), ones(1,p)], [zeros(1,q), linspace(0,1,ncp-q+2), ones(1,q)]};
%! Lx  = 1; Ly  = 1;
%! [cntl(1,:,:), cntl(2,:,:)] = meshgrid(linspace(0, Lx, ncp+1), linspace(0, Ly, mcp+1) );
%! cntl(4,:,:) = 1:numel(cntl(1,:,:));
%! nrb = nrbmak(cntl, knots);
%! u = rand (1, 30); v = rand (1, 10);
%! N = nrbnumbasisfun ({u, v}, nrb);
%! assert (all(all(N>0)), true)
%! assert (all(all(N<=(ncp+1)*(mcp+1))), true)
%! assert (max(max(N)),(ncp+1)*(mcp+1))
%! assert (min(min(N)),1)