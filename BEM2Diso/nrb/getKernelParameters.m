function [ jacob, normals, r,relDist ] = getKernelParameters(crv, dcrv,collocCoords, xi )
 % calculate the parameters we need to evaluate the kernels
[fieldPt,dxydxi] = nrbdeval(crv, dcrv, xi) ;

% dxydxi=dN*elcoords;             % the geometry derivatives
jacob=norm(dxydxi);             % jacobian
normals=1/jacob * [ dxydxi(2) -dxydxi(1) ];

% fieldPt=N*elcoords;
relDist=fieldPt(1:2)'-collocCoords;
r=norm(relDist);
% dr=1/r * relDist;
% drdn=dr(1:2)*normals';

end

