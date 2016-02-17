clear all
init;

for i = 1:mc.nrun  
    
    tic;   
     %% update zHH
    [z_HH_Individuals,z_HH] = UpdateHouseholdIndicator(...
        orig.n,orig.n_individuals,orig.maxd,orig.p, hyper.K,hyper.L,...
    para.phi,orig.dataT,orig.SS,para.w,para.pi,para.HHdata_all,para.lambda);

    %% update zmember
    z_Individuals = UpdateIndividualIndicator(para.phi,orig.dataT,para.w,...
        hyper.K,hyper.L,orig.p,orig.maxd,z_HH,orig.HHserial);
 
    %% generate impossible house hold
    [z_Individual_extra,z_HH_extra,para.IndividualData_all,para.HHdata_all,para.hh_size_new] = ...
    GetImpossibleHouseholds(para.lambda,para.w,para.phi,orig.d,orig.p,hyper.L, orig.ACS_count,para.pi,...
    orig.origdata,orig.HHdataorig);
    
    %% combine data and indicators
    para.z_HH_all = [z_HH;z_HH_extra];
    para.z_HH_Individuals_all = [z_HH_Individuals;z_Individual_extra(:,1)];  
    para.z_Individual_all = [z_Individuals;z_Individual_extra(:,2)];
    clear z_HH z_HH_extra z_Individuals z_HH_Individuals z_Individual_extra
    
    %% update phi
    [para.phi,para.phicountcluster,para.kcount] = UpdatePhi(...
        hyper.K,hyper.L,orig.p,orig.d,orig.maxd,...
        para.IndividualData_all,para.z_Individual_all,para.z_HH_all,para.z_HH_Individuals_all);
    
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
   
