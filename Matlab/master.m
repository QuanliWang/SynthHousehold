clear all
init;

for i = 1:nrun  
    
    tic;   
    n_household_all = length(HHdata1_all); %total number of household 
    
    %% update zHH
    [z_HH_all, z_HH_Individuals_all,newphi] = UpdateHouseholdIndicator(n,n_individuals,maxd,p, K,L,...
    phi,dataT,SS,w,pi,HHdata1_all,lambda1,HHdata2_all,lambda2,ImpossibleIndividuals,z_HH_extra);

    %% update zmember
    z_Individual_all = UpdateIndividualIndicator(n,n_individuals,...
    newphi,dataT,w,K,L,p,maxd,z_HH_all,HHserial,ImpossibleIndividuals);

    %% update phi
    [phi,phicountcluster,kcount,phicount] = UpdatePhi(data_full_all,K,L,p,d,maxd,...
        z_HH_Individuals_all,z_Individual_all,z_HH_all);
    
    %% update lambda
    [lambda1,lambda2] = UpdateLambda(K,dHH1,dHH2,z_HH_all,HHdata1_all,HHdata2_all);
    
    %% update pi
    [pi,u] = UpdatePi(K,kcount,alpha);
        
    %% update w
    [w,v] = UpdateW(K,L,phicountcluster,beta);
     
    %% update alpha
    alpha = UpdateAlpha(aa,ab,K,u);
    
    %% update alpha
    beta = UpdateBeta(K,L,v,ba,bb);
    
    %% generate impossible house hold
    [hh_size_new,ImpossibleIndividuals,data_full_all,z_HH_extra,HHdata1_all,HHdata2_all] = ...
    GetImpossibleHouseholds(lambda1,lambda2,w,phi,d,p,L, ACS_count,pi,...
    origdata,HHdataorig);
   
    postsave;
    
end
   
