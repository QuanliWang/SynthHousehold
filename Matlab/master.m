clear all
init;

for i = 1:nrun  
    
    tic;   
    n_household_all = size(HHdata_all,1); %total number of household 
    
    %% update zHH
    [z_HH_all, z_HH_Individuals_all,newphi] = UpdateHouseholdIndicator(...
        orig.n,orig.n_individuals,orig.maxd,orig.p, hyper.K,hyper.L,...
    phi,orig.dataT,orig.SS,w,pi,HHdata_all,lambda1,lambda2,ImpossibleIndividuals,z_HH_extra);

    %% update zmember
    z_Individual_all = UpdateIndividualIndicator(orig.n,orig.n_individuals,...
    newphi,orig.dataT,w,hyper.K,hyper.L,orig.p,orig.maxd,z_HH_all,orig.HHserial,ImpossibleIndividuals);

    %% update phi
    [phi,phicountcluster,kcount] = UpdatePhi(data_full_all,hyper.K,hyper.L,orig.p,orig.d,orig.maxd,...
        z_HH_Individuals_all,z_Individual_all,z_HH_all);
    
    %% update lambda
    [lambda1,lambda2] = UpdateLambda(hyper.dHH,hyper.K,z_HH_all,HHdata_all);
    
    %% update pi
    [pi,u] = UpdatePi(alpha,kcount);
        
    %% update w
    [w,v] = UpdateW(beta,phicountcluster);
     
    %% update alpha
    alpha = UpdateAlpha(hyper.aa,hyper.ab,u);
    
    %% update alpha
    beta = UpdateBeta(hyper.ba,hyper.bb,v);
    
    %% generate impossible house hold
    [hh_size_new,ImpossibleIndividuals,data_full_all,z_HH_extra,HHdata_all] = ...
    GetImpossibleHouseholds(lambda1,lambda2,w,phi,orig.d,orig.p,hyper.L, orig.ACS_count,pi,...
    orig.origdata,orig.HHdataorig);
   
    postsave;
    
end
   
