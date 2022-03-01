function p = Ftest2tail(F,dfn,dfd)

% function p = Ftest2tail(F,dfn,dfd)
% 
% calculates the two-tailed probability of finding a larger value of F, P(>F);
% dfn and dfd are the numbers of degrees of freedom
% for the numerator and denominator, respectively

p = 2.0*betainc(dfd./(dfd+dfn*F),0.5*dfd,0.5*dfn);
if any(p>1.0)
    p2 = 2.0 - p;
    p = min([p p2], [], 2);
end
