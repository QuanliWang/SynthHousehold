function beta = UpdateBeta(K,L,v,ba,bb)
    vnew = v(:,1:L-1);
    beta = gamrnd(ba + K*(L-1), 1/(bb - sum(sum(log(1-vnew)))));
end

