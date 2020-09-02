% Given a spacing, coil diameters, and integer of complete turns; return
% Exp[L] and arbitrarily precise trace Width.
function [L, w] = currentSheet(s, n, Dout, Din)
    
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
   
   u = pi*(4e-7);     % Permeability of free space.
   
   L = 0.5*c1*u*n*n*Davg*(log(c2/rho)+c3*rho+c4*rho^2)*1e-6;

endfunction