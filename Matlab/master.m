clear all
init;

for i = 1:mc.nrun  
    
    tic;   
     %% update zHH
    [z_HH_Individuals,z_HH] = UpdateHouseholdIndicator(...
        orig.n,orig.n_individuals,orig.maxd,orig.p, orig.dataT,orig.SS,...
        hyper.K,hyper.L,...
        para.phi,para.w,para.pi,para.HHdata_all(1:orig.n,:),para.lambda);

    %% update zIndividual
    z_Individuals = UpdateIndividualIndicator(orig.p,orig.maxd,orig.dataT,...
        orig.HHserial, hyper.K,hyper.L,z_HH,para.phi,para.w);
 
    %% generate impossible house hold
    [z_Individual_extra,z_HH_extra,IndividualData_extra,HHdata_extra,para.hh_size_new] = ...
    GetImpossibleHouseholds(orig.d,orig.ACS_count,...
                            para.lambda,para.w,para.phi, para.pi);
    
    %% combine data and indicators
    para.z_HH_all = [z_HH;z_HH_extra];
    para.HHdata_all = orig.HHdataorig; para.HHdata_all(:,2) = para.HHdata_all(:,2) - 1;
    para.HHdata_all = [para.HHdata_all;HHdata_extra];
    para.IndividualData_all = [orig.origdata(:,1:8);IndividualData_extra(:,1:8)];
    %column 1 for K groups and column 2 for L groups
    para.z_Individual_all = [[z_HH_Individuals;z_Individual_extra(:,1)]...
                               [z_Individuals;z_Individual_extra(:,2)]];
                           
    clear z_HH z_HH_extra z_Individuals z_HH_Individuals z_Individual_extra IndividualData_extra HHdata_extra
    
    %% update phi
    para.phi = UpdatePhi(hyper.K,hyper.L,orig.p,orig.d,orig.maxd,...
                para.IndividualData_all,para.z_Individual_all);
    
     %% update lambda
    [para.lambda] = UpdateLambda(hyper.dHH,hyper.K,para.z_HH_all,para.HHdata_all);
    
     %% update pi
    [para.pi,para.u] = UpdatePi(para.alpha,para.z_HH_all,hyper.K);
        
    %% update w
    [para.w,para.v] = UpdateW(para.beta,para.z_Individual_all, hyper.K, hyper.L);
    
     %% update alpha
    para.alpha = UpdateAlpha(hyper.aa,hyper.ab,para.u);
    
    %% update beta
    para.beta = UpdateBeta(hyper.ba,hyper.bb,para.v);
    
    postsave;
    
end
   
