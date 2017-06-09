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

function [B, id] = nrbbasisfun (points, nrb)

% NRBBASISFUN: Basis functions for NURBS
%
% Calling Sequence:
% 
%   B      = nrbbasisfun (u, crv)
%   B      = nrbbasisfun ({u, v}, srf)
%   [B, N] = nrbbasisfun ({u, v}, srf)
%   [B, N] = nrbbasisfun (p, srf)
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
%      B - Basis functions 
%          size(B)=[numel(u),(p+1)] for curves
%          or [numel(u)*numel(v), (p+1)*(q+1)] for surfaces
%
%      N - Indices of the basis functions that are nonvanishing at each
%          point. size(N) == size(B)
%   

  if (~iscell(nrb.knots))    %% NURBS curve
    
    [B, id] = nrb_crv_basisfun__ (points, nrb);
    
  else                       %% NURBS surface

    if (iscell(points))
      [v, u] = meshgrid(points{2}, points{1});
      p(1,:,:) = u;
      p(2,:,:) = v;
      p = reshape(p, 2, []);
    else
      p = points;
    end
    
    [B, id] = nrb_srf_basisfun__ (p, nrb); 

  end
end

  function [B, nbfu] = nrb_crv_basisfun__ (points, nrb);
    n    = size (nrb.coefs, 2) -1;
    p    = nrb.order -1;
    u    = points;
    U    = nrb.knots;
    w    = nrb.coefs(4,:);
    
    spu  =  findspan (n, p, u, U);
    nbfu =  numbasisfun (spu, u, p, U);
    
    N     = w(nbfu+1) .* basisfun (spu, u, p, U);
    B     = bsxfun (@(x,y) x./y, N, sum (N,2));

  end
  

%!test
%! U = [0 0 0 0 1 1 1 1];
%! x = [0 1/3 2/3 1] ;
%! y = [0 0 0 0];
%! w = rand(1,4);
%! nrb = nrbmak ([x;y;y;w], U);
%! u = linspace(0, 1, 30);
%! B = nrbbasisfun (u, nrb);
%! xplot = sum(bsxfun(@(x,y) x.*y, B, x),2);
%!
%! yy = y; yy(1) = 1;
%! nrb2 = nrbmak ([x.*w;yy;y;w], U); 
%! #figure, plot(xplot, B(:,1), nrbeval(nrb2, u)(1,:).', w(1)*nrbeval(nrb2, u)(2,:).')
%! assert(B(:,1), w(1)*nrbeval(nrb2, u)(2,:).', 1e-6)
%! 
%! yy = y; yy(2) = 1;
%! nrb2 = nrbmak ([x.*w;yy;y;w], U);
%! #figure, plot(xplot, B(:,2), nrbeval(nrb2, u)(1,:).', w(2)*nrbeval(nrb2, u)(2,:).')
%! assert(B(:,2), w(2)*nrbeval(nrb2, u)(2,:).', 1e-6)
%!
%! yy = y; yy(3) = 1;
%! nrb2 = nrbmak ([x.*w;yy;y;w], U);
%! #figure, plot(xplot, B(:,3), nrbeval(nrb2, u)(1,:).', w(3)*nrbeval(nrb2, u)(2,:).')
%! assert(B(:,3), w(3)*nrbeval(nrb2, u)(2,:).', 1e-6)
%!
%! yy = y; yy(4) = 1;
%! nrb2 = nrbmak ([x.*w;yy;y;w], U);
%! #figure, plot(xplot, B(:,4), nrbeval(nrb2, u)(1,:).', w(4)*nrbeval(nrb2, u)(2,:).')
%! assert(B(:,4), w(4)*nrbeval(nrb2, u)(2,:).', 1e-6)

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
%! [B, N] = nrbbasisfun ({u, v}, nrb);
%! assert (sum(B, 2), ones(300, 1), 1e-6)
%! assert (all(all(B<=1)), true)
%! assert (all(all(B>=0)), true)
%! assert (all(all(N>0)), true)
%! assert (all(all(N<=(ncp+1)*(mcp+1))), true)
%! assert (max(max(N)),(ncp+1)*(mcp+1))
%! assert (min(min(N)),1)