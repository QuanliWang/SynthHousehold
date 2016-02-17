clear all
init;

for i = 1:mc.nrun  
    
    tic;   
    
    %% generate impossible house hold
    [para.hh_size_new,z_Individual_extra,para.data_full_all,z_HH_extra,para.HHdata_all] = ...
    GetImpossibleHouseholds(para.lambda,para.w,para.phi,orig.d,orig.p,hyper.L, orig.ACS_count,para.pi,...
    orig.origdata,orig.HHdataorig);
    
    %% update zHH
    [para.z_HH_all, para.z_HH_Individuals_all] = UpdateHouseholdIndicator(...
        orig.n,orig.n_individuals,orig.maxd,orig.p, hyper.K,hyper.L,...
    para.phi,orig.dataT,orig.SS,para.w,para.pi,para.HHdata_all,para.lambda,...
    z_Individual_extra,z_HH_extra);
    
    %% update zmember
    para.z_Individual_all = UpdateIndividualIndicator(orig.n,orig.n_individuals,...
    para.phi,orig.dataT,para.w,hyper.K,hyper.L,orig.p,orig.maxd,para.z_HH_all,...
    orig.HHserial,z_Individual_extra);

    %para.n_household_all = size(para.HHdata_all,1); %total number of household
    
    %% update phi
    [para.phi,para.phicountcluster,para.kcount] = UpdatePhi(para.data_full_all,...
        hyper.K,hyper.L,orig.p,orig.d,orig.maxd,...
        para.z_HH_Individuals_all,para.z_Individual_all,para.z_HH_all);
    
     %% update lambda
    [para.lambda] = UpdateLambda(hyper.dHH,hyper.K,para.z_HH_all,para.HHdata_all);
    
     %% update pi
    [para.pi,para.u] = UpdatePi(para.alpha,para.kcount);
        
    %% update w
    [para.w,para.v] = UpdateW(para.beta,para.phicountcluster);
    
     %% update alpha
    para.alpha = UpdateAlpha(hyper.aa,hyper.ab,para.u);
    
    %% update beta
    para.beta = UpdateBeta(hyper.ba,hyper.bb,para.v);
    
    postsave;
    
end
   
