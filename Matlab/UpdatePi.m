function [pi,u] = UpdatePi(alpha,K,kcount)
    %tic;
    u = zeros(K,1);
    pi = zeros(K,1);
    for k = 1:K-1
        u(k) = betarnd(1 + kcount(k),alpha + sum(kcount(k+1:K)));
        if u(k)>1-1e-5
            u(k)=1-1e-5;
        end  
    end
    u(K) = 1;

    pi(1:K)=u.*cumprod([1;1-u(1:K-1)]);
end

