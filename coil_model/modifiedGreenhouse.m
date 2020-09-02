% Given a spacing, coil diameters, and integer of complete turns; return
% Exp[L] and arbitrarily precise trace Width.
function [L, w] = modifiedGreenhouse(s, n, Dout, Din, t)
    
   % Useful definitions to carry through:
   delta = (Dout - Din)/2;
   Davg  = (Dout + Din)/2;
   rho   = (Dout - Din) / (Dout + Din);
   
   % Calculate the width from the given (integer) number of turns:
   w = (delta - s*(n-1))/n;
   
   % Using Modified Wheeler approximation from Mohan's 1999 paper.
   c1 = 1.27;       % Table 2: c1 Square.
   c2 = 2.07;       % Table 2: c2 Square.
   c3 = 0.18;       % Table 2: c3 Square.
   c4 = 0.13;       % Table 2: c4 Square.
   
   u = pi*(4e-7)*1e-6;     % Permeability of free space.
   mu = u;
   len = 4*n*Dout - 4*n*w - (2*n-1)^2*(s+w) + w    ;
   if n>1
     d_ = (w+s) * (3*n - 2*n - 1) * (n+1) / 3 / (n-1);
     L = u/2/pi*len*(log(len/n/(w+t))-0.2-0.47*n+(n-1)*(log(sqrt(1+(len/4/n/d_)^2)+len/4/n/d_)-sqrt(1+(4*n*d_/len)^2)+4*n*d_/len)) 
   else
     L = u/2/pi*len*(log(len/n/(w+t))-0.2-0.47*n);
   end
%{
   L = 0.5*c1*u*n*n*Davg*(log(c2/rho)+c3*rho+c4*rho^2);


   
   s_ = delta / (2*n-1);
   w_ = (delta - s_*(n-1))/n;
   len_ = 4*n*Dout - 4*n*w_ - (2*n-1)^2*(s_+w_) + w_    ;
   L_w_ = u/2/pi*len_*log(len_/n/(w_+t));
   
   L_w = u/2/pi*len*log(len/n/(w+t));
   
   L = L + L_w - L_w_;
%}
   
endfunction