function H = Hz(s, n, Dout, Din, n_q, y_q)

H = 0;
w = ((Dout - Din)/2 - s*(n-1))/n;

for i = 1 : n
  D(i) = Dout - 2*(i-1)*s-2*(i-0.5)*w;
end

for i = 1 : n
  r = (s+w)*(n_q-i)-y_q;
  b = D(i)/2;
  a = -D(i)/2;
  H += H_segment(w,r,a,b);
  H += H_segment(w,D(i)-r,a,b);
  H += H_segment(w,-a,-r,D(i)-r);
  H += H_segment(w,b,r-D(i),r);
end


function H_ = H_segment(w_,r_,a_,b_)
  if abs(r_)<w_/2
    H_=0;
  else
    H_ = 1 / (4*pi*r_) * (b/sqrt(b^2+r_^2)-a/sqrt(a^2+r_^2));
  endif
%  H_ = r;
%  H_ = integral2(@(x,y) 1 / (4*pi*w_) ./ (x.^2+(y-r_).^2).^1.5 .* (y-r_), a_, b_, -w_/2, w_/2,"AbsTol",1e-3,"RelTol",1e-3);
%  H_ = 0.5*(log(sqrt(b_^2+(r_+w_)^2)-b)-log(sqrt(b_^2+(r_+w_)^2)+b)+log(sqrt(a_^2+(r_+w_)^2)+a)-log(sqrt(a_^2+(r_+w_)^2)-a)+log(sqrt(b_^2+(r_-w_)^2)+b)-log(sqrt(b_^2+(r_-w_)^2)-b)+log(sqrt(a_^2+(r_-w_)^2)-a)-log(sqrt(a_^2+(r_-w_)^2)+a))
endfunction

endfunction