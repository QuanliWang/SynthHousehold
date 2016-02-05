function data_to_check = GenerateData2Check(hh_size,lambda1, lambda2, w, ...
    phi,pi, d, p, number_of_generation, cum_number_of_generation,L)
    data_to_check = zeros(10000,8*hh_size+1+hh_size);
    for m=1:10000
        pi_size = pi.*lambda2(:,hh_size-1);
        pi_size_renorm = pi_size./sum(pi_size);
        hhindexh = randomsample(pi_size_renorm,rand);  
        data_to_check(m,hh_size * 8 + 1) = hhindexh;

        syn = zeros(hh_size,p+1);
        for hh = 1:hh_size
            memberindexhh = randomsample(w(hhindexh,:),rand);
            data_to_check(m,hh_size * 8 +1 + hh) = memberindexhh;
            % generating individual level data
            for j = 1:p
                phimj = reshape(phi((hhindexh-1)*L+memberindexhh,j,1:d(j)),1,d(j));
                syn(hh,j) = randomsample(phimj,rand);
            end
        end
        % generating household level data
        lambda1mj = lambda1(hhindexh,:);
        syn(:,p+1) = randomsample(lambda1mj,rand);
        syn_sorted = sortrows(syn,5);

        for hh = 1:hh_size
            ibase = (hh-1) * 8;
            data_to_check(m,ibase + 1) = m + (cum_number_of_generation+number_of_generation)*10000;
            data_to_check(m,ibase + 2) = hh;
            data_to_check(m, ibase + (3:8)) = syn_sorted(hh,:);    
        end
        
    end        
end

