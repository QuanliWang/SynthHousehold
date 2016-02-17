 %%    
if (mod(i,mc.thin) == 0 && i > mc.burn) 
   output.piout((i-mc.burn)/mc.thin,:) = pi';
   output.wout((i-mc.burn)/mc.thin,:,:) = w;
   output.newphiout((i-mc.burn)/mc.thin,:,:) = newphi;
   output.lambda1out((i-mc.burn)/mc.thin,:,:) = lambda1;
   output.lambda2out((i-mc.burn)/mc.thin,:,:) = lambda2;
end
i
n_household_all
toc
output.elapsed_time(i) = toc;
output.n_sout(i) = size(data_full_all,1);
output.nout(i) = n_household_all;
output.extrasize(i,:) = hh_size_new(2:4);
output.z_HH_save(i,1:orig.n_individuals) = z_HH_Individuals_all(1:orig.n_individuals);
output.z_member_save(i,1:orig.n_individuals) = z_Individual_all(1:orig.n_individuals);
output.alphaout(i) = alpha;
output.betaout(i) = beta;