clear all
init;

for i = 1:nrun  
    
    tic;   
    n_household_all = size(HHdata_all,1); %total number of household 
    
    %% update zHH
    [z_HH_all, z_HH_Individuals_all,newphi] = UpdateHouseholdIndicator(...
        orig.n,orig.n_individuals,orig.maxd,orig.p, K,L,...
    phi,orig.dataT,orig.SS,w,pi,HHdata_all,lambda1,lambda2,ImpossibleIndividuals,z_HH_extra);

    %% update zmember
    z_Individual_all = UpdateIndividualIndicator(orig.n,orig.n_individuals,...
    newphi,orig.dataT,w,K,L,orig.p,orig.maxd,z_HH_all,orig.HHserial,ImpossibleIndividuals);

    %% update phi
    [phi,phicountcluster,kcount] = UpdatePhi(data_full_all,K,L,orig.p,orig.d,orig.maxd,...
        z_HH_Individuals_all,z_Individual_all,z_HH_all);
    
    %% update lambda
    [lambda1,lambda2] = UpdateLambda(dHH,K,z_HH_all,HHdata_all);
    
    %% update pi
    [pi,u] = UpdatePi(alpha,kcount);
        
    %% update w
    [w,v] = UpdateW(beta,phicountcluster);
     
    %% update alpha
    alpha = UpdateAlpha(aa,ab,K,u);
    
    %% update alpha
    beta = UpdateBeta(ba,bb,v);
    
    %% generate impossible house hold
    [hh_size_new,ImpossibleIndividuals,data_full_all,z_HH_extra,HHdata_all] = ...
    GetImpossibleHouseholds(lambda1,lambda2,w,phi,orig.d,orig.p,L, orig.ACS_count,pi,...
    orig.origdata,orig.HHdataorig);
   
    postsave;
    
end
   
