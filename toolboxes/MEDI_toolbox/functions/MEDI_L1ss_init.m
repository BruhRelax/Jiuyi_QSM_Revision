function [chi_p_init,chi_n_init,R2p, alpha, beta] = MEDI_L1ss_init(Mask,CF,R2s,QSM,delta_TE)
       
        alpha = double(1.21/3*CF/1e6)*2*pi;
        beta = double(1/3*CF/1e6)*2*pi;
        N=sum(Mask(:)>0);
        RR=real(double(R2s(Mask>0))/100);
        QQ=real(double(QSM(Mask>0)));
        ub = [4*ones(N,1);zeros(N,1)];
        lb = [zeros(N,1);-4*ones(N,1)];
        x0=zeros(2*N,1);
        options = optimoptions('lsqnonlin',...
            'Algorithm','trust-region-reflective',...
            'SpecifyObjectiveGradient', true,...
            'Display', 'off',...
            'JacobianMultiplyFcn',@MEDI_L1ss_init_jmfcn);
        C2=@(x)MEDI_L1ss_init_fcn(x,alpha,beta,RR,QQ);
        x = lsqnonlin(C2,x0,lb,ub,options);
        zpp = zeros(size(Mask)); 
        zpp(Mask>0)=x(1:N);
        znp = zeros(size(Mask)); 
        znp(Mask>0)=x(N+1:2*N);
        

        chi_p_init = zpp*(2*pi*delta_TE*CF)/1e6;
        chi_n_init = znp*(2*pi*delta_TE*CF)/1e6;
        R2p = R2s;
end
