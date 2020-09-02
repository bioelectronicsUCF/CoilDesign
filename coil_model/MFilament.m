function M = MFilament(rxsp, rxn, rxod, rxid, txsp, txn, txod, txid, z)
  rxw = ((rxod-rxid)/2-rxsp*(rxn-1))/rxn;
  txw = ((txod-txid)/2-txsp*(txn-1))/txn;
  mu = 4 * pi * 1e-7 * 1e-6;
  tx_w_seg_n=1;
  rx_w_seg_n=1;
  M = 0;
  for tx_seg = 1 : txn
    for rx_seg = 1 : rxn
      for tx_w_seg = 1 : tx_w_seg_n
        for rx_w_seg = 1 : rx_w_seg_n
          a = txod/2 - (tx_seg-1)*(txw+txsp) - (tx_w_seg-1)/tx_w_seg_n*txw;
          b = rxod/2 - (rx_seg-1)*(rxw+rxsp) - (rx_w_seg-1)/rx_w_seg_n*rxw;
          M_=0;
          M_ = M_+sqrt(2*(a+b)^2+z^2);
          M_ = M_+sqrt(2*(a-b)^2+z^2);
          M_ = M_-2*sqrt(2*a^2+2*b^2+z^2);
          M_ = M_-(a+b)*atanh((a+b)/sqrt(2*(a+b)^2+z^2));
          M_ = M_-(a-b)*atanh((a-b)/sqrt(2*(a-b)^2+z^2));
          M_ = M_+(a+b)*atanh((a+b)/sqrt(2*a^2+2*b^2+z^2));
          M_ = M_+(a-b)*atanh((a-b)/sqrt(2*a^2+2*b^2+z^2));
          M = M + 1/(tx_w_seg_n*rx_w_seg_n)*2*mu/pi*M_;
        endfor
      endfor
    endfor
  endfor
    
  M = abs(M);

endfunction