function [w,v] = UpdateW(beta,z_Individual_all, K, L)
    
    phicountcluster = groupcount(z_Individual_all(1,:),z_Individual_all(2,:),K,L);
    cum = fliplr(cumsum(fliplr(phicountcluster),2));
    
    v = betarnd(1 + phicountcluster(:,1:(L-1)),beta + cum(:,2:L));
    v(v > 1-1e-5) = 1 - 1-1e-5;
    v = [v  ones(K,1)];
    
    w = v.* cumprod([ones(K,1) 1 - v(:,1:(L-1))],2);

end

