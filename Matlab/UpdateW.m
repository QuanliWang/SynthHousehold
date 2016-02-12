function [w,v] = UpdateW(K,L,phicountcluster,beta)
    %tic;
    v = zeros(K,L);
    w=zeros(K,L);
    for k = 1:K
        vk = v(k,:);
        for l = 1:L-1
            vk(l) = betarnd(1 + phicountcluster(k,l), beta + sum(phicountcluster(k,l+1:L)));
            if vk(l)>1-1e-5
                vk(l)=1-1e-5;
            end
            vk(L) = 1;
        end
        v(k,:) = vk;
        w(k,1:L) = vk.*cumprod([1,1-vk(1:L-1)]);
    end
    %toc
    disp('w updated');
end

