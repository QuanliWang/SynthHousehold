init;

for i = 1:nrun  
    
    UpdateParameters;
    
    [hh_size_new,trueextras,data_full,z_HH_new,HHdata1_all,HHdata2_all] = ...
    GetImpossibleHouseholds(lambda1,lambda2,w,phi,d,p,L, ACS_count,pi,...
    origdata,HHdataorig);
   
    postsave;
    
end
   
