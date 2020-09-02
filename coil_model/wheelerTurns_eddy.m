% Given a spacing, coil diameters, and integer of complete turns; return
% Exp[L] and arbitrarily precise trace Width.
function [L, w] = wheelerTurns_eddy(s, n, Dout, Din, t, f, Sheet_R_coil)
    
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
   
   len = 4*n*Dout - 4*n*w - (2*n-1)^2*(s+w) + w    ;
   L_w = u/2/pi*len*log(len/n/(w+t));
   
   L = L + L_w - L_w_;   
   
   Hc = n * 1 / (Dout-Din) * log ( ((Dout^2+t^2)^0.5+Dout) / ((Din^2+t^2)^0.5+Din) );

   Hin = Hc * exp(Din/Dout*(0.4+0.15*log(Dout/2/t)));
   Hout = -Hc * (0.1+0.08*log(Dout/2/t)) * exp(Din/Dout*(1+1/8*log(Dout/2/t)));
  
   H = Hout - (linspace(0,(n-1)/n,n) + 1/2/n) * (Hout-Hin);
   Hin_ = Hout - linspace(1/n,1,n) * (Hout-Hin);
   Hout_ = Hout - linspace(0,(n-1)/n,n) * (Hout-Hin);
   delH = (Hin-Hout)/n;
   delH_ = linspace(1,1,n) * delH;
   
   omega = 2 * pi * f;
   rho = Sheet_R_coil*t;
   mu = 4 * pi * 1e-7 * 1e-6;
   delta = (rho/pi/f/mu)^0.5;

   [I_eddy_edge, I_eddy_center] = Ieddy(w,delta,H,Sheet_R_coil,omega,n,t);

   for i = 1:n
     M_edge(i) = MFilament2(s,n,Dout,Din,0,1,Dout-2*(i-1)*(w+s),Dout-2*(i-1)*(w+s)-2*delta,0)-MFilament2(s,n,Dout,Din,0,1,Dout-2*(i-1)*(w+s)-2*w+2*delta,Dout-2*(i-1)*(w+s)-2*w,0);
     M_center(i) = MFilament2(s,n,Dout,Din,0,1,Dout-2*(i-1)*(w+s)-2*delta,Dout-2*(i-1)*(w+s)-w,0)-MFilament2(s,n,Dout,Din,0,1,Dout-2*(i-1)*(w+s)-w,Dout-2*(i-1)*(w+s)-2*w+2*delta,0);
     L = L - M_edge(i) * real(I_eddy_edge(i)) - M_center(i) * real(I_eddy_center(i));
   end  
     
end