 %%    
if (mod(i,mc.thin) == 0 && i > mc.burn) 
   output.piout((i-mc.burn)/mc.thin,:) = para.pi';
   output.wout((i-mc.burn)/mc.thin,:,:) = para.w;
   output.newphiout((i-mc.burn)/mc.thin,:,:) = para.phi;
   output.lambda1out((i-mc.burn)/mc.thin,:,:) = para.lambda{1};
   output.lambda2out((i-mc.burn)/mc.thin,:,:) = para.lambda{2};
end
i
total_household = size(para.HHdata_all,1)
toc
output.elapsed_time(i) = toc;
%output.n_sout(i) = size(para.IndividualData_all,1);
output.nout(i) = total_household;
output.extrasize(i,:) = para.hh_size_new;
output.z_HH_save(i,:) = para.z_Individual_all(1:orig.n_individuals,1);
output.z_member_save(i,:) =para. z_Individual_all(1:orig.n_individuals,2);
output.alphaout(i) = para.alpha;
output.betaout(i) = para.beta;