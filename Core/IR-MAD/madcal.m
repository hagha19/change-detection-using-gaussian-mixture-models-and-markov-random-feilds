function [Xcal theta ntrain ntest] = madcal(Xuncal, Yref, nochangeProb, nochangeProbPerc)

%
% *** TEST VERSION ***
%
% [Xcal theta ntrain ntest] = madcal(Xuncal, Yref, nochangeProb, nochangeProbPerc);
%
% MADCAL reads output from IRMADENVI and performs (orthogonal) regression
%   on no-change observations for normalization
%
% MADCAL appends a report to file 'regreport.txt'
% If 'regreport.txt' does not exist it is created
%
% Input
% Xuncal            - uncalibrated data (3-D)
% Yref              - reference data (3-D)
% nochangeProb      - probability of no-change (as output w from program irmadenvi)
% nochangeProbPerc  - percentile (optional, defaults to 95%)
%
% Output
% Xcal              - calibrated data
% theta             - calibration coefficients for bands
% ntrain            - number of training pixels
% ntest             - number of test pixels

% (c) Copyright 2008-2010
% Allan Aasbjerg Nielsen
% aa@space.dtu.dk, www.imm.dtu.dk/~aa
% 20 Sep 2010

if nargin < 3, error('madcal: too few input arguments'); end
if nargin < 4, nochangeProbPerc = 0.95; end
if nargin > 4, error('madcal: too many input arguments'); end

% Proportion of observations set aside for testing
testproportion = 0.3333333;
orthoreg = 1; % For OLS use orthoreg=0, for orthogonal regression use orthoreg=1
seed = 4711007; % for initialization of rand

readfile1 = 0; % Xuncal is 3-D Matlab matrix
readfile2 = 0; % Yref is 3-D Matlab matrix
if ischar(Xuncal), readfile1 = 1; end % Xuncal is filename
if ischar(Yref), readfile2 = 1; end % Yref is filename

if readfile1 == 1
    [Xuncal,p] = freadenvit(Xuncal); % from file
    nr = p(2); nc = p(1); nv = p(3);
else
    if ndims(Xuncal)~=3, error('input1 must be 3-D'); end
    [nr, nc, nv] = size(Xuncal);
end
if readfile2 == 1
    [Yref,p] = freadenvit(Yref); % from file
    nry = p(2); ncy = p(1); nvy = p(3);
else
    if ndims(Yref)~=3, error('input1 must be 3-D'); end
    [nry, ncy, nvy] = size(Yref);
end

if (nr ~= nry) | (nc ~= ncy) | (nv ~= nvy)
    error('madcal: Xuncal and Yref do not match')
end
if ndims(nochangeProb)~=2, error('nochangeProb must be 2-D'); end
[nrp, ncp] = size(nochangeProb);
if (nr ~= nrp) | (nc ~= ncp)
    error('madcal: Xuncal/Yref and nochangeProb do not match')
end

Xuncal = reshape(Xuncal, nr*nc, nv);
Yref = reshape(Yref, nr*nc, nv);
nochangeProb = reshape(nochangeProb, nr*nc, 1);

% hard threshold nochangeProb image to find no-change pixels
xxx = Xuncal(nochangeProb > nochangeProbPerc, :);
yyy = Yref(nochangeProb > nochangeProbPerc, :);
clear Yref nry ncy nvy nrp ncp

[nobs nvar] = size(xxx);

% set aside some observations for training and some for testing
rand('state',seed);
idx = rand(nobs, 1) > testproportion;
% training data
x  = xxx(idx, :);
y  = yyy(idx, :);
nobs = size(x, 1);
ntrain = nobs;
% test data
xt = xxx(~idx, :);
yt = yyy(~idx, :);
nobst = size(xt, 1);
ntest = nobst;

band = 1:nvar;

%% correlations
R = corrcoef([x y]);
R = diag(R(nvar+1:end,1:nvar))';
% Rt = corrcoef([xt yt]);
% Rt = diag(Rt(nvar+1:end,1:nvar))';

v1 = ones(nobs,1);
v1t = ones(nobst,1);

thetao = [];
thetah = [];
% thetaw = [];
rmse = [];
stderr = [];
corr = [];
t = [];
pt = [];
yth = [];
% rmset = [];
% chisq = [];
% pchisq = [];

diary 'regreport.txt';

%% Print report on regression analyses
% disp(' ');
disp('Number of training observations:');
disp(int2str(nobs));
disp(' ');

disp('Orthogonal regressions for normalization (traning observations only)');
disp(['Channel Intercept Std.err. t p Slope Std.err. t p Correlation RMSE'])

for i=1:nvar
    xx = x(:,i);
    yy = y(:,i);
    %% Orthogonal regression start
    [c,b] = eig(cov([xx yy]));
    aux = c(1,1)/c(2,1);
    thetao  = [thetao [aux*mean(xx)+mean(yy); -aux]];
    alpha = thetao(1, end);
    beta  = thetao(2, end);
    
    xmean = mean(xx);
    ymean = mean(yy);
    covmat = cov([xx yy],1); % divide by nobs instead of (nobs-1)
    sxx = covmat(1,1);
    syy = covmat(2,2);
    sxy = covmat(1,2);
%     % estimates from Kendall and Stuart (1979), vol. 2, p. 405 and p. 402
%     beta = (syy-sxx+sqrt((sxx-syy)^2+4*sxy^2))/(2*sxy);
%     alpha = ymean-beta*xmean;
    
    % sigma2 from Bilbo (1989) p. 56
    %%sigma2 = 1/4*(sxx+syy-sqrt((sxx+syy)^2-4*(sxx*syy-sxy^2)));
    %sigma2 = 1/4*(sxx+syy-sqrt((sxx-syy)^2+4*sxy^2));
    % sigma2 from Kendall and Stuart (1979), vol. 2, p. 411
    sigma2 = 1/(2*(1+beta^2))*(syy-2*beta*sxy+beta^2*sxx);
    sigma2 = 2*nobs/(nobs-2)*sigma2;
    % dispersion of [alpha beta]' from Patefield (1977)
    tauhat = sigma2*beta/((1+beta^2)*sxy);
    aux = (1+beta^2)*sigma2*beta/(nobs*sxy);
    aux1 = xmean*(1+tauhat);
    dispmat = [xmean*aux1+sxy/beta -aux1; -aux1 1+tauhat];
    dispmat = aux.*dispmat;
    
%     X = [v1 xx];
%     yh = X*thetao(:,end);
    yh = alpha + beta*xx;
    eh = yy-yh;
    df = ntrain-2;
    mse = eh'*eh/df;
    rmse = [rmse sqrt(mse)];
    stderr = [stderr sqrt(diag(dispmat))];
    t = thetao(:, i) ./ stderr(:,i);
    pt = ttest(t,df);
%     pt = betainc(df./(df+t.^2),0.5*df,0.5);
    
    %% Orthogonal regression stop

    % %% Schott et al. regression start
    %     aux = sqrt(var(yy)/var(xx));
    %     thetaw = [thetaw [mean(yy)-aux*mean(xx); aux]];
    % %% Schott et al. regression stop

    %% OLS
    X = [v1 xx];
    thetah = [thetah X\yy];
    if orthoreg ~= 1
        yh = X*thetah(:,end);
        eh = yy-yh;
        df = nobs-size(X,2);
        mse = eh'*eh/df;
        rmse = [rmse sqrt(mse)];
        dispersion = rmse(end)*inv(X'*X);
        stderr = [stderr sqrt(diag(dispersion))];
        corr = [corr diag(1./stderr(:,end))*dispersion*diag(1./stderr(:,end))];
        t = [t thetah(:,end)./stderr(:,end)];
        pt = [pt betainc(df./(df+t(:,end).^2),0.5*df,0.5)];
    end

    %% Print report on regression analyses
    disp([int2str(i) '   ' ...
            num2str(alpha)  '   ' ...
            num2str(stderr(1,end)) '   ' num2str(t(1)) '   ' num2str(pt(1)) '   ' ...
            num2str(beta)  '   ' ...
            num2str(stderr(2,end)) '   ' num2str(t(2)) '   ' num2str(pt(2)) '   ' ...
            num2str(R(i)) '   ' num2str(rmse(i))])
%            num2str(R(i)) ' & ' num2str(rmse(i)) ' & ' ...
%            num2str(Rt(i)) ' & ' num2str(rmset(i)) '\\'])

    figure, plot(xx,yy,'k.');
    hold on
    xrange = [0 max(xx)];
    yrange = [0 max(yy)];
    range = [xrange yrange];
    plot(xrange, thetah(1,end)+xrange*thetah(2,end),'k:'); % OLS
    plot(xrange, thetao(1,end)+xrange*thetao(2,end),'k-'); % orthogonal regression
%     plot(xrange, thetaw(1,end)+xrange*thetaw(2,end),'k-.'); % Schott et al. regression
    %axis(range);
    axis equal;
    title(['Band ' num2str(i)]);
    xlabel('uncalibrated');
    ylabel('reference');
    hold off
    print('-deps2',['fig',int2str(i),'.eps'])

    figure, plot(yh, yy, 'k.');
    hold on
    xrange = [0 max(yh)];
    yrange = [0 max(yy)];
    range = [xrange yrange];
    plot(xrange, xrange,'k-'); %y=x
    %axis(range);
    axis equal;
    title(['Band ' num2str(i) ', training data']);
    ylabel('reference');
    xlabel('calibrated');
    hold off
    print('-deps2',['calibfigtrain',int2str(i),'.eps'])
    
    %% examine test observations after transformation
    if (orthoreg == 1)
        yth = [yth [v1t xt(:,i)]*thetao(:,end)];
    else
        yth = [yth [v1t xt(:,i)]*thetah(:,end)];
    end
% %       yth = [yth [v1t xt(:,i)]*thetaw(:,end)]; % Schott et al. regression
%     eth = yt(:,i)-yth(:,end);
%     dft = size(yth(:,end),1);
%     mset = eth'*eth/dft;
%     rmset = [rmset sqrt(mset)];
    if orthoreg
        Xcal(:, i) = thetao(1, i) + Xuncal(:, i)*thetao(2, end);
    else
        Xcal(:, i) = thetah(1, i) + Xuncal(:, i)*thetah(2, end);
    end    
end
Xcal = reshape(Xcal, nr, nc, nv);

if orthoreg
    theta = thetao;
else
    theta = thetah;
end

%% Test of test observations for MAD radiometric calibration

disp(' ')
disp('Number of test observations:');
disp(int2str(nobst));

m0 = mean(xt)';
m1 = mean(yt)';
m2 = mean(yth)';
S0 =cov(xt);
S1 = cov(yt);
S2 = cov(yth);

%% Individual scalar paired t-test for equal means
% dif = yt-yth;
dif = yth-yt;
mdif = mean(dif)';
Smdif = cov(dif)./nobst;
df3 = nobst-1;
tpair = mdif./sqrt(diag(Smdif));
ptpair = ttest(tpair,df3);
%ptpair = betainc(df3./(df3+tpair.^2),0.5*df3,0.5);% P(>|tpair|)
disp(' ');
disp('Paired t-tests for equal means, we want p-values > 0.05 (the closer to 1, the better)');
disp('------------------------------');
disp('channel:');
disp(int2str(band));
disp('Uncorrected:');
disp(num2str(m0'));
disp('Normalized:')
disp(num2str(m2'));
disp('Reference:')
disp(num2str(m1'));
disp('Difference:');
disp(num2str(mdif'));
disp('t:');
disp(num2str(tpair'));
% disp('%P{>|t|}:');
disp('p:');
disp(num2str(ptpair'));

%% Individual scalar F-tests for equal variances (assuming unknown means)
var0 = diag(S0);
var1 = diag(S1);
var2 = diag(S2);
F = var2./var1;
% FF = F;
% for i=1:nvar
%     if F(i)<1;
%         FF(i)=1/F(i);
%     end
% end
% pF = 2*Ftest(FF,nobst-1,nobst-1);
pF2 = Ftest2tail(F,nobst-1,nobst-1);
disp(' ');
disp('F-tests for equal variances, we want p-values > 0.05 (the closer to 1, the better)');
disp('---------------------------');
disp('channel:');
disp(int2str(band));
disp('Uncorrected:');
disp(num2str(var0'));
disp('Normalized:')
disp(num2str(var2'));
disp('Reference:');
disp(num2str(var1'));
disp('F:');
disp(num2str(F'));
% % disp('%If F>1 P{<=1/F}+P{>F}; if F<1 P{<=F}+P{>1/F}:')
% disp('p:');
% disp(num2str(pF'));
disp('p:');
disp(num2str(pF2'));

disp(' ');
disp(' ');

diary off

% return

for i=1:nvar
    figure, plot(yth(:, i),yt(:, i),'k.');
    hold on
    xrange = [0 max(yth(:,i))];
    yrange = [0 max(yt(:,i))];
    range = [xrange yrange];
    plot(xrange, xrange,'k-'); %y=x
    %axis(range);
    axis equal;
    title(['Band ' num2str(i) ', test data']);
    ylabel('reference');
    xlabel('calibrated');
    hold off
    print('-deps2',['calibfigtest',int2str(i),'.eps'])
end

% return

figure, errorbar(theta(1,:), stderr(1,:), '.-')
if orthoreg
    title('Intercepts, orthogonal regression')
else
    title('Intercepts, ordinary regression')
end

figure, errorbar(theta(2,:), stderr(2,:), '.-')
if orthoreg
    title('Slopes, orthogonal regression')
else
    title('Slopes, ordinary regression')
end

figure, plot(ptpair, '.-')
ylim([0 1]);
title('Paired t-tests for equal means, p-values')

figure, plot(pF2, '.-')
ylim([0 1]);
title('F-tests for equal variances, p-values')
