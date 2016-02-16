 %%    
if (mod(i,thin) == 0 && i > burn) 
   piout((i-burn)/thin,:) = pi';
   wout((i-burn)/thin,:,:) = w;
   newphiout((i-burn)/thin,:,:) = newphi;
   lambda1out((i-burn)/thin,:,:) = lambda1;
   lambda2out((i-burn)/thin,:,:) = lambda2;
end
i
n_household_all
toc
elapsed_time(i) = toc;
n_sout(i) = size(data_full_all,1);
nout(i) = n_household_all;
extrasize(i,:) = hh_size_new(2:4);
z_HH_save(i,1:orig.n_individuals) = z_HH_Individuals_all(1:orig.n_individuals);
z_member_save(i,1:orig.n_individuals) = z_Individual_all(1:orig.n_individuals);
alphaout(i) = alpha;
betaout(i) = beta;