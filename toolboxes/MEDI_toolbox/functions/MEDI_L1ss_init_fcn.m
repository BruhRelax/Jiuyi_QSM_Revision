function [Y,Jinfo] = MEDI_L1ss_init_fcn(x,alpha,beta,RR,QQ)
rep=@(y)repmat(y, [1 size(x,2)]);
% https://www.mathworks.com/support/search.html/answers/1899800-when-using-jacobianmultiplyfcn-in-lsqnonlin-why-is-jinfo-required-to-be-numeric-how-can-i-pass-mor.html
s.alpha=alpha;
s.beta=beta;
Jinfo=JMInfo(s);
N=size(x,1)/2;
Y=[alpha*x(1:N,:)/100-beta*x(N+1:2*N,:)/100-rep(RR);x(1:N,:)+x(N+1:2*N,:)-rep(QQ)];
end

