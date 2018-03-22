function calcula_arco(x1,y1,x2,y2,xc,yc)
  # Function to compute the tet1 angle between the line from point (x1,y1) to (xc,yc) and the
  # horizontal direction and the the tet2 angle between the line from point (x2,y2) to (xc,yc)
  # and the horizontal direction


  dx1 = x1 - xc; dy1 = y1 - yc;
  dx2 = x2 - xc; dy2 = y2 - yc;

  # Computation of tet1
  if dy1 == 0				# The point 1 and the center have the same y coordinate
    if x1 > xc
      tet1 = 0;
    else  # (x1 < xc)
      tet1 = pi;
    end;
  elseif dx1 == 0				# The point 1 and the center have the same x coordinate
    if y1 > yc
      tet1 = pi/2;
    else  # (y1 < yc)
      tet1 = -pi/2;
    end;
  else  # (dx1~=0 e dy1~=0)
    tet1 = atan(dy1/dx1);
    if dx1<0 && tet1<0
      tet1 = pi + tet1;
    elseif dx1 < 0 && tet1>0
      tet1 = -pi + tet1;
    end;
  end;

  # Computation of tet2
  if dy2 == 0				# The point 2 and the center have the same y coordinate
    if x2 > xc
      tet2 = 0;
    else  # (x2 < xc)
      tet2 = pi;
    end;
  elseif dx2 == 0				# The point 2 and the center have the same x coordinate
    if y2 > yc
      tet2 = pi/2;
    else  # (y2 < yc)
      tet2 = -pi/2;
    end;
  else  # (dx2~=0 e dy2~=0)
    tet2 = atan(dy2/dx2);
    if dx2<0 && tet2<0
      tet2 = pi + tet2;
    elseif dx2 < 0 && tet2>0
      tet2 = -pi + tet2;
    end;
  end;
  return tet1, tet2
end
