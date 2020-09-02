function output=p_deriv(var, point)
  args = point(:);
  dx = 1e-6;
  y1 = hypergeom([args(1),args(2)],args(3),args(4));
  args(var) = args(var) + dx;
  y2 = hypergeom([args(1),args(2)],args(3),args(4));
  output = (y1-y2)/dx;
end