% Given a spacing, coil diameters, and integer of complete turns; return
% Exp[L] and arbitrarily precise trace Width.
function [L, w] = wheelerTurns(s, n, Dout, Din, t)
    
   % Useful definitions to carry through:
   delta = (Dout - Din)/2;
   Davg  = (Dout + Din)/2;
   rho   = (Dout - Din) / (Dout + Din);
   
   % Calculate the width from the given (integer) number of turns:
   w = (delta - s*(n-1))/n;
   
   % Using Modified Wheeler approximation from Mohan's 1999 paper.
   K1 = 2.34;       % Table 1: K1 Square.
   K2 = 2.75;       % Table 1: K2 Square.
   u = pi*(4e-7)*1e-6;     % Permeability of free space.
   
   L = K1*u*n*n*Davg/(1+K2*rho);
 

   s_ = delta / (2*n-1);
   w_ = (delta - s_*(n-1))/n;
   len_ = 4*n*Dout - 4*n*w_ - (2*n-1)^2*(s_+w_) + w_    ;
   L_w_ = u/2/pi*len_*log(len_/n/w_);
   
   len = 4*n*Dout - 4*n*w - (2*n-1)^2*(s+w) + w  ;
   L_w = u/2/pi*len*log(len/n/(w+t));
   
   L = L + abs(L_w) - abs(L_w_);
end