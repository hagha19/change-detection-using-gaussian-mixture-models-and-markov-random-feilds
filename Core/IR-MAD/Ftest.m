function p = Ftest(F,dfn,dfd)

% function p = Ftest(F,dfn,dfd)
% 
% calculates the probability of finding a larger value of F, P(>F);
% dfn and dfd are the numbers of degrees of freedom
% for the numerator and denominator, respectively

p = betainc(dfd./(dfd+dfn*F),0.5*dfd,0.5*dfn); % P(>F)