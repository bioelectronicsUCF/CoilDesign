function [I_eddy_edge, I_eddy_center] = Ieddy_iterartion(w,delta,H,Sheet_R_coil,omega,n,t,s)
  
 [I_eddy_edge, I_eddy_center] = Ieddy(w,delta,H,Sheet_R_coil,omega,n,t);
  for i = 1:n
    for j = 1:n
      if (abs(i - j)<1)
      elseif (abs(i-j)==1)
        H(i) = H(i) - I_eddy_edge(j)/2/pi * (((log(1/(s+w+delta/2))-log(1/(s+delta/2)))/w)-(log(1/(s+w*2-delta/2))-log(1/(s+w-delta/2)))/w);
        H(i) = H(i) - I_eddy_center(j)/2/pi * (((log(1/(s+w+w/4))-log(1/(s+w/4)))/w)-(log(1/(s+w*2-w/4))-log(1/(s+w-w/4)))/w);
      endif
    end
  end
 [I_eddy_edge, I_eddy_center] = Ieddy(w,delta,H,Sheet_R_coil,omega,n,t);

endfunction
