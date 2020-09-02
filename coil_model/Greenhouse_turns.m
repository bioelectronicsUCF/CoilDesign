function L_turns = Greenhouse_turns(s, n, od, id, t);

  w = ((od-id)/2 - s*(n-1))/n;
  u = pi*(4e-7)*1e-6;     % Permeability of free space.

  for i = 1 : n
    len = (4*(od-w)-8*(i-1)*(w+s));
    L_self_turns(i) = u/2/pi*len*( log(len/(w+t)) + 0.2 );
  end
  
  M_minus_turns = linspace(0,0,n);
  M_plus_turns = linspace(0,0,n);
  for i = 1 : n
    for j = 1 : n
      M_minus_turns(i) = M_minus_turns(i) + 2 * M_seg_(i,j,"-");
      M_plus_turns(i) = M_plus_turns(i) + 2 * M_seg_(i,j,"+"); 
    endfor
  endfor
  
  L_turns = L_self_turns - M_minus_turns + M_plus_turns;
    
  function Q_ = Q_cal(_len,_gmd)
    if (_len==0) || (_gmd==0)
      Q_=0;
    else
      Q_ = log( _len/_gmd + sqrt( 1+(_len/_gmd)^2) ) + _gmd/_len - sqrt( 1+(_gmd/_len)^2 );
    end
  endfunction
  
  function M_ij_seg_ = M_seg_(i_,j_,sign_)
 
    if (sign_ == "-")||(sign_ == "+")
      
      len_seg_i_ = od-w-2*(i_-1)*(w+s);
      len_seg_j_ = od-w-2*(j_-1)*(w+s);
      len_avg_ = (len_seg_i_+len_seg_j_)/2;
      len_delta_ = abs(len_seg_i_-len_seg_j_)/2;
      if sign_ == "-"
        gmd_ = len_avg_;
      elseif sign_ == "+"
        gmd_ = w+s;
      endif
      
      Q_avg_ = Q_cal(len_avg_,gmd_);
      M_avg_ = u/pi * len_avg_ * Q_avg_;
      Q_delta_ = Q_cal(len_delta_,gmd_);
      M_delta_ = u/pi * len_delta_ * Q_delta_;
      M_ij_seg_ = M_avg_ - M_delta_;

    else
      M_ij_seg_ = 0;
    endif    
    
  endfunction
  
endfunction