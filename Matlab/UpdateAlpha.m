function alpha = UpdateAlpha(aa,ab,u)
    K = length(u);
    alpha = gamrnd(aa + K - 1, 1/(ab - sum(log(1-u(1:K-1)))));
end
