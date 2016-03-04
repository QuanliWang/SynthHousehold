%%%%---define output files---%%
output.alphaout = zeros(mc.nrun,1); 
output.betaout = zeros(mc.nrun,1);
output.piout = zeros(mc.eff.sam,hyper.K); 
output.wout = zeros(mc.eff.sam,hyper.K,hyper.L);
output.nout = zeros(mc.nrun,1);

output.extrasize = zeros(mc.nrun,3);

output.z_HH_save = zeros(mc.nrun,orig.n_individuals);
output.z_member_save = zeros(mc.nrun,orig.n_individuals);
output.elapsed_time = zeros(1,mc.nrun);

output.newphiout = zeros(mc.eff.sam,orig.maxd*orig.p,hyper.K*hyper.L);
output.lambda1out = zeros(mc.eff.sam,hyper.K,hyper.dHH(1));
output.lambda2out = zeros(mc.eff.sam,hyper.K,hyper.dHH(2));