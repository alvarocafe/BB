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
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

%% -*- texinfo -*-
%% @deftypefn {Function File} {[@var{dRdx}, @var{dRdy}]=} nrbbasisfungradient (@{@var{dzdu}, @var{dzdv}@}, @{@var{dxdu}, @var{dydu}, @var{dxdv}, @var{dydv}@})
%% @deftypefnx {Function File} {[@var{dRdx}]=} nrbbasisfungradient (@var{dzdu}, @var{dxdu})
%% NRBBASISFUNGRADIENT Compute the gradient of the basis functions of a NURBS surface at the
%% specified parametric points.
%%
%% INPUT:
%% @itemize @minus
%% @item @var{dzdu}, @var{dzdv}: basis function derivatives with respect
%% to parameters @code{u} and @code{v}
%% @item @var{dxdu}, @var{dydu}, @var{dxdv},
%% @var{dydv}: NURBS geometry map derivatives
%% @end itemize
%%
%% OUTPUT:
%% @itemize @minus
%% @item @var{dRdx} derivative of the basis functions with respect to
%% the @code{x} direction
%% @item @var{dRdy} derivative of the basis functions with respect to
%% the @code{y} direction
%% @end itemize
%%
%% @seealso{nrbbasisfunder,nrbbasisfun,nrbderiv}
%% @end deftypefn

%% Author: Carlo de Falco <carlo@guglielmo.local>
%% Created: 2009-04-29

function [varargout]  = nrbbasisfungradient (dz, N, nrb)

  if (iscell(dz))

    dzdu=dz{1}; 
    dzdv=dz{2}; 

    X = squeeze(nrb.coefs(1,:,:)./nrb.coefs(4,:,:));
    Y = squeeze(nrb.coefs(2,:,:)./nrb.coefs(4,:,:));

    nfun = size(N,2);
    tmp = sum(X(N).*dzdu, 2);
    dxdu= tmp(:, ones(nfun,1)); 
    tmp = sum(Y(N).*dzdu, 2);
    dydu= tmp(:, ones(nfun,1));
    tmp = sum(X(N).*dzdv, 2);
    dxdv= tmp(:, ones(nfun,1));
    tmp = sum(Y(N).*dzdv, 2);
    dydv= tmp(:, ones(nfun,1));
    
    detjac = dxdu.*dydv - dydu.*dxdv;
    
    varargout{1} = ( dydv .* dzdu - dydu .*dzdv)./detjac;
    varargout{2} = (-dxdv .* dzdu + dxdu .*dzdv)./detjac;
    %%keyboard
  elseif (~iscell(dz) && (nargout==1))

    varargout{1} = dzdu ./ dxdu;

  else
    
    print_usage();
    
  end

end

%!test
%! 
%! p = 2;
%! q = 3;
%! mcp = 2; ncp = 3;
%! knots = {[zeros(1,p), linspace(0,1,mcp-p+2), ones(1,p)], [zeros(1,q), linspace(0,1,ncp-q+2), ones(1,q)]};
%! Lx  = 1; Ly  = 1;
%! [cntl(2,:,:), cntl(1,:,:)] = meshgrid(linspace(0, Ly, ncp+1), linspace(0, Lx, mcp+1) );
%! cntl(4,:,:) = 1:numel(cntl(1,:,:));
%! nrb = nrbmak(cntl, knots);
%! 
%! [u(2,:,:), u(1,:,:)] = meshgrid(rand (1, 20), rand (1, 20));
%! 
%! 
%! uv (1,:) = u(1,:,:)(:)';
%! uv (2,:) = u(2,:,:)(:)';
%! 
%! [dzdu, dzdv, connect] =  nrbbasisfunder (uv, nrb);
%! nd        = nrbderiv(nrb);
%! [ndp, dp] = nrbdeval(nrb, nd, uv);
%! 
%! dxdu = repmat (reshape (dp{1}(1,:,:), [], 1), 1, columns(dzdu));
%! dydu = repmat (reshape (dp{1}(2,:,:), [], 1), 1, columns(dzdu));
%! dxdv = repmat (reshape (dp{2}(1,:,:), [], 1), 1, columns(dzdu));
%! dydv = repmat (reshape (dp{2}(2,:,:), [], 1), 1, columns(dzdu));
%! 
%! [dzdx, dzdy]= nrbbasisfungradient ({dzdu, dzdv}, {dxdu, dydu, dxdv, dydv});
%! assert(norm(sum(dzdx, 2), inf), 0, 1e-10)
%! assert(norm(sum(dzdy, 2), inf), 0, 1e-10)