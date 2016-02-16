function [SS,origdata00] = PrepareData(file,Household)
    size234alldata = load(file);
    newdata = size234alldata.newdata;

    HHindex = unique(newdata(:,1)); 
    SSfull = hist(newdata(:,1),HHindex); 

    randindex = randsample(length(SSfull),Household);
    SS = SSfull(randindex);
    
    origdata00 = [];
    for i=1:Household
        index_i = find(newdata(:,1)==randindex(i));
        base_vec = ones(length(index_i),1);
        currentHH = [base_vec * i newdata(index_i,2:end) base_vec * SS(i)];
        origdata00 = [origdata00;currentHH];
    end
    
end

