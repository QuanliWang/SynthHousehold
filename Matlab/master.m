clear all
init;

%synindex = [4960 4970 4980 4990 5000 9960 9970 9980 9990 10000];
synindex = [10 15 20 25];
synData = cell(length(synindex),1);
for i = 1:mc.nrun  
    
    tic;   
     %% update zHH
     [z_HH,z_HH_Individuals] = samplezHHwithHHnewv1_2HHvar(para.phi,orig.dataT,...
         para.w,para.pi,orig.SS,para.HHdata_all(:,1:orig.n)',para.lambda{1},para.lambda{2}, rand(orig.n,1));

     %% update zIndividual
     z_Individuals = samplezmemberv1(para.phi,orig.dataT,para.w,z_HH,orig.HHserial,rand(size(orig.dataT,2),1));
    

        %% generate impossible house hold
        [z_Individual_extra,z_HH_extra,IndividualData_extra,HHdata_extra,para.hh_size_new,current_synData] = ...
        GetImpossibleHouseholds(orig.d,orig.ACS_count,para.lambda,para.w,...
        para.phi,para.pi,hyper.howmany,orig.n, ismember(i,synindex));
        if ismember(i,synindex) 
            synData{find(synindex ==i)} = current_synData(1:8,:);
        end
        %% combine data and indicators
        para.z_HH_all = [z_HH z_HH_extra];
        para.HHdata_all = orig.HHdataorigT; para.HHdata_all(2,:) = para.HHdata_all(2,:) - 1;
        para.HHdata_all = [para.HHdata_all HHdata_extra];
        para.IndividualData_all = [orig.origdata(:,1:8)' IndividualData_extra];
        
        %row 1 for K groups and row 2 for L groups
        para.z_Individual_all = [[z_HH_Individuals;z_Individuals] z_Individual_extra];
       
        clear z_HH z_HH_extra z_Individuals z_HH_Individuals z_Individual_extra IndividualData_extra HHdata_extra

        %% update phi
        para.phi = UpdatePhi(hyper.K,hyper.L,orig.p,orig.d,orig.maxd,...
                    para.IndividualData_all,para.z_Individual_all);
        %% update w
        [para.w,para.v] = UpdateW(para.beta,para.z_Individual_all, hyper.K, hyper.L);
  
    
     %% update lambda
    [para.lambda] = UpdateLambda(hyper.dHH,hyper.K,para.z_HH_all,para.HHdata_all);
    
     %% update pi
    [para.pi,para.u] = UpdatePi(para.alpha,para.z_HH_all,hyper.K);
        
    
    
     %% update alpha
    para.alpha = UpdateAlpha(hyper.aa,hyper.ab,para.u);
    
    %% update beta
    para.beta = UpdateBeta(hyper.ba,hyper.bb,para.v);
    
    postsave;
    
%    if mod(i,1000) == 0
%        save -v7.3 last;
%     end
end
   

