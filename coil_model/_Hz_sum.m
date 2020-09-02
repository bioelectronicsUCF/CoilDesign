function _H = _Hz_sum(s, n, Dout, Din, n_q, y_q, n_k)
  w = ( (Dout-Din)/2 - (n-1)*s )/n;
  
  _H = 0;
  _n_k = n_k;
  _x = Din/2 + w/2 + (s+w)*(n-n_q) + y_q;
  for _k = 1 : _n_k
    _w = w / _n_k;
    _s = s + w - _w;
    _n = n;
    _Dout = Dout - 2*(_k-1)*_w;
    _Din = _Dout - 2*(_n-1)*_s - 2*_n*_w;
    if _x>_Dout/2
      _n_q = 1;
      _y_q = _x-_Dout/2+_w/2;
    elseif _x<_Din/2+_w+_s/2
      _n_q = _n;
      _y_q = _x-_Din/2-_w/2;
    else
      _n_q = _n - floor((_x-_Din/2+_s/2)/(_s+_w));
      _y_q = _x - (_Din/2+_w/2+(_n-_n_q)*(_s+_w));
    endif
    _H+= Hz(_s, _n, _Dout, _Din, _n_q, _y_q)/_n_k;
  endfor
endfunction