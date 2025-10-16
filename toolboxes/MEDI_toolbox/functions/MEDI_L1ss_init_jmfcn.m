function w = MEDI_L1ss_init_jmfcn(Jinfo,Y,flag)
%JMFCN Summary of this function goes here
%   Detailed explanation goes here

% https://www.mathworks.com/support/search.html/answers/1899800-when-using-jacobianmultiplyfcn-in-lsqnonlin-why-is-jinfo-required-to-be-numeric-how-can-i-pass-mor.html
alpha=Jinfo.s.alpha;
beta=Jinfo.s.beta;

N=size(Y,1)/2;
if flag > 0
    w = Jpositive(Y);
elseif flag < 0
    w = Jnegative(Y);
else
    w = Jnegative(Jpositive(Y));
end
    function w = Jpositive(x)
        w=[alpha*x(1:N,:)/100-beta*x(N+1:2*N,:)/100;x(1:N,:)+x(N+1:2*N,:)];
    end
    function w = Jnegative(x)
        w=[alpha*x(1:N,:)/100+x(N+1:2*N,:);-beta*x(1:N,:)/100+x(N+1:2*N,:)];
    end
end

