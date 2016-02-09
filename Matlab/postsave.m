 %%    
    disp('synthesis updated');
    
       %syndata{i} = syndatao;
       if (mod(i,thin) == 0 && i > burn) 
           piout((i-burn)/thin,:) = pi';
           wout((i-burn)/thin,:,:) = w;
           newphiout((i-burn)/thin,:,:) = newphi;
           lambda1out((i-burn)/thin,:,:) = lambda1;
           lambda2out((i-burn)/thin,:,:) = lambda2;
       end
       i
       n_new
       toc
       elapsed_time(i) = toc;
       n_sout(i) = n_s_new;
       nout(i) = n_new;
       size2extrasize(i) = hh_size_new(2);
       size3extrasize(i) = hh_size_new(3);
       size4extrasize(i) = hh_size_new(4);
       z_HH_save(i,1:n_s) = z_HH_all;
       z_member_save(i,1:n_s) = z_member;
       alphaout(i) = alpha;
       betaout(i) = beta;
       disp('saving done'); 