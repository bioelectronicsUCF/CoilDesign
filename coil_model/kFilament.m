function k = kFilament(rxsp, rxn, rxod, rxid, txsp, txn, txod, txid, z, rxt, txt, f, rxrsheet, txrsheet)

  L_rx = wheelerTurns_eddy(rxsp,rxn,rxod,rxid,rxt,f,rxrsheet);
  L_tx = wheelerTurns_eddy(txsp,txn,txod,txid,txt,f,txrsheet);
  
  rho = (4/pi)^(1+min([rxod,txod])/max([rxod,txod]));
  mu = 4 * pi * 1e-7 * 1e-6;  
  txw = ((txod-txid)/2+txsp)/(txn-1)-txsp;
  rxw = ((rxod-rxid)/2+rxsp)/(rxn-1)-rxsp;
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
  
  k = M / sqrt(L_rx*L_tx);

end