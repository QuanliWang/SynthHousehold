 %%    
if (mod(i,thin) == 0 && i > burn) 
   output.piout((i-burn)/thin,:) = pi';
   output.wout((i-burn)/thin,:,:) = w;
   output.newphiout((i-burn)/thin,:,:) = newphi;
   output.lambda1out((i-burn)/thin,:,:) = lambda1;
   output.lambda2out((i-burn)/thin,:,:) = lambda2;
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