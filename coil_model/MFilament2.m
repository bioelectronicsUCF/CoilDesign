function M = MFilament2(rxsp, rxn, rxod, rxid, txsp, txn, txod, txid, z)

  txw = ((txod-txid)/2+txsp)/txn-txsp;
  rxw = ((rxod-rxid)/2+rxsp)/rxn-rxsp;

  rho = (4/pi)^(1+rxod/txod);
  mu = 4 * pi * 1e-7 * 1e-6;  %H/um
  
  M = 0;
  for tx_seg = 1 : txn
    for rx_seg = 1 : rxn
      a = txod/2 - (tx_seg-1)*(txw+txsp) - txw/2;
      b = rxod/2 - (rx_seg-1)*(rxw+rxsp) - rxw/2;
      gamma = 2 * a * b / (a^2 + b^2 + z^2);
      
      M_ = mu * pi * a^2 * b^2 / (2 * (a^2 + b^2 + z^2) ^ 1.5) * (1 + 15/32 * gamma^2 + 315/1024 *gamma^4);
      M = M + rho * M_;
    end
  end

end