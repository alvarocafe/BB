function nurbs = nrbmak(coefs,knots)
%
% Function Name:
%
%   nrbmak - Construct the NURBS structure given the control points
%            and the knots.
%
% Calling Sequence:
%
%   nurbs   = nrbmak(cntrl,knots);
%
% Parameters:
%
%   cntrl       : Control points, these can be either Cartesian or
% 		homogeneous coordinates.
%
% 		For a curve the control points are represented by a
% 		matrix of size (dim,nu) and for a surface a multidimensional
% 		array of size (dim,nu,nv). Where nu is number of points along
% 		the parametric U direction, and nv the number of points
%               along the V direction. Dim is the dimension valid options
% 		are
% 		2 .... (x,y)        2D Cartesian coordinates
% 		3 .... (x,y,z)      3D Cartesian coordinates
% 		4 .... (wx,wy,wz,w) 4D homogeneous coordinates
%
%   knots	: Non-decreasing knot sequence spanning the interval
%               [0.0,1.0]. It's assumed that the curves and surfaces
%               are clamped to the start and end control points by knot
%               multiplicities equal to the spline order.
%               For curve knots form a vector and for a surface the knot
%               are stored by two vectors for U and V in a cell structure.
%               {uknots vknots}
%
%   nurbs 	: Data structure for representing a NURBS curve.
%
% NURBS Structure:
%
%   Both curves and surfaces are represented by a structure that is
%   compatible with the Spline Toolbox from Mathworks
%
% 	nurbs.form   .... Type name 'B-NURBS'
% 	nurbs.dim    .... Dimension of the control points
% 	nurbs.number .... Number of Control points
%       nurbs.coefs  .... Control Points
%       nurbs.order  .... Order of the spline
%       nurbs.knots  .... Knot sequence
%    nurbs.uKnots    .... unique knots
%   Note: the control points are always converted and stored within the
%   NURBS structure as 4D homogeneous coordinates. A curve is always stored
%   along the U direction, and the vknots element is an empty matrix. For
%   a surface the spline degree is a vector [du,dv] containing the degree
%   along the U and V directions respectively.
%


nurbs.form   = 'B-NURBS';
nurbs.dim    = 4;
np = size(coefs);
dim = np(1);


% constructing a curve
nurbs.number = np(2);
if (dim < 4)
    nurbs.coefs = repmat([0.0 0.0 0.0 1.0]',[1 np(2)]);
    nurbs.coefs(1:dim,:) = coefs;
else
    nurbs.coefs = coefs;
end
nurbs.order = size(knots,2)-np(2);
knots = sort(knots);
nurbs.knots = (knots-knots(1))/(knots(end)-knots(1));
nurbs.uKnots=unique(knots);

ne=length(nurbs.uKnots)-1;   % number of elements
elRange=zeros(ne,2);        % initialise matrices
elConn=zeros(ne,nurbs.order);
elKnotIndices=zeros(ne,2);

element=1;
previousKnotVal=0;
for i=1:length(knots)
    currentKnotVal=knots(i);
    if knots(i)~=previousKnotVal
        elRange(element,:)=[previousKnotVal currentKnotVal];
        elKnotIndices(element,:)=[i-1 i];
        element=element+1;
    end
    previousKnotVal=currentKnotVal;
end
nurbs.elRange=elRange;
nurbs.elKnotIndices=elKnotIndices;

numRepeatedKnots=0;
for e=1:ne
    indices=(elKnotIndices(e,1)-nurbs.order+2):elKnotIndices(e,1);
    previousKnotVals=knots(indices);
    currentKnotVals=ones(1,nurbs.order-1)*knots(elKnotIndices(e,1));
    if isequal(previousKnotVals,currentKnotVals) && length(nonzeros(previousKnotVals))>1;
        numRepeatedKnots=numRepeatedKnots+1;
    end
    elConn(e,:)=(elKnotIndices(e,1)-nurbs.order+1):elKnotIndices(e,1);
end
elConn(end,end)=1;  % the last point is equal to the first point
nurbs.ELEM=elConn;

numBasisFns=length(knots)-nurbs.order;
collocPts=zeros(1,numBasisFns);
for i=1:numBasisFns
    collocPts(i)=sum(knots((i+1):(i+(nurbs.order-1))))/(nurbs.order-1);
end
nurbs.collocPts=collocPts;
