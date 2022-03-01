function p = ttest(t,df)

% function p = ttest(t,df)
% 
% calculates the probability of finding a larger value of the absolute value of t,
% P(>|t|);
% df is the number of degrees of freedom

p=betainc(df./(df+t.^2),0.5*df,0.5);