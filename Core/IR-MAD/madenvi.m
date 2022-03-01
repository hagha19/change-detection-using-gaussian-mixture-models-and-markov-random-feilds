function [mads,rho,v1,v2,s11,s22,s12,prob,chi2] = madenvi(X,Y)

% [mads,rho,v1,v2,s11,s22,s12] = madenvi(X,Y)
%
% MAD - multivariate alteration detection
%
% Input
% X      - file name of multivariate band sequential input image
%          number one or image already read into Matlab
% Y      - file name of multivariate band sequential input image
%          number two or image already read into Matlab
% 
% Output
% the MAD variates (mads),
% the canonical correlations (rho),
% the eigenvectors (v1 and v2) normed to give CVs unit variance, and
% the relevant (variance-)covariance matrices, and
% if desired, the probability of no-change and the approximate chi-squared quantity. 

% (c) Copyright 2005-2010
% Allan Aasbjerg Nielsen, Ph.D., M.Sc.
% aa@space.dtu.dk, www.imm.dtu.dk/~aa
% 19 Sep 2010

if nargin<2, error('Not enough input arguments.'); end
if nargin>2, error('Too many input arguments.'); end
% if nargin==3
%     fnameo = varargin{1};
%     if ~ischar(fnameo), error('fnameo should be a char string'); end
% end

readfile1 = 0; % input1 is 3-D matrix
readfile2 = 0; % input2 is 3-D matrix
if ischar(X), readfile1 = 1; end % input1 is filename
if ischar(Y), readfile2 = 1; end % input2 is filename

if readfile1==1
    [X,p] = freadenvit(X);
    nrows = p(2); ncols = p(1); nvar1 = p(3);
else
    if ndims(X)~=3, error('input1 must be 3-D'); end
    [nrows,ncols,nvar1] = size(X);
end
if readfile2==1
    [Y,p] = freadenvit(Y);
    nrow2 = p(2); ncol2 = p(1); nvar2 = p(3);
else
    if ndims(Y)~=3, error('input2 must be 3-D'); end
    [nrow2,ncol2,nvar2] = size(Y);
end
X = reshape(X,nrows*ncols,nvar1);
Y = reshape(Y,nrow2*ncol2,nvar2);

if nvar2>nvar1
    error('input with highest number of variables must be first set');
end
if ~(ncols==ncol2 & nrows==nrow2)
    error('data in fname1 and fname2 do not match');
end

N = nrows*ncols;
covxy = cov([X Y]);
s11 = covxy(1:nvar1,1:nvar1);
s22 = covxy(nvar1+1:end,nvar1+1:end);
s12 = covxy(1:nvar1,nvar1+1:end);
s21 = s12';

if nvar2==nvar1 % solve smallest eigenproblem

%[v1,d1] = eigen2(s12*(s22^(-1))*s21,s11);
invs22 = inv(s22);
[v1,d1] = eigen2(s12*invs22*s21,s11);
rho = diag(sqrt(d1))'; % lowest canonical correlation first

% scale v1 to give CVs with unit variance
aux1 = v1'*s11*v1; % dispersion of CVs
aux2 = 1./sqrt(diag(aux1));
aux3 = repmat(aux2',nvar1,1);
v1 = v1.*aux3; % now dispersion is unit matrix
%v1'*s11*v1

% sum of correlations between X and CV(X) positive
invstderr = diag(1./std(X));
%invstderrcv = diag(1./sqrt(diag(v1'*s11*v1)));
%sgn = diag(sign(sum(invstderr*s11*v1*invstderrcv)));
sgn = diag(sign(sum(invstderr*s11*v1)));
v1 = v1*sgn;
%figure; bar(sum(invstderr*s11*v1))

%%[v2,d2] = eigen2(s21*(s11^(-1))*s12,s22);
%[v2,d2] = eigen2(s21*inv(s11)*s12,s22);
%v2 = v2*diag(sign(diag(v1'*s12*v2)));
v2 = invs22*s21*v1; %./repmat(rho,nvar1,1); % scaling doesn't matter

% scale v2 to give CVs with unit variance
aux1 = v2'*s22*v2; % dispersion of CVs
aux2 = 1./sqrt(diag(aux1));
aux3 = repmat(aux2',nvar2,1);
v2 = v2.*aux3; % now dispersion is unit matrix
%v2'*s22*v2
%invstderr = diag(1./std(Y)); figure; bar(sum(invstderr*s22*v2))

mads = (X-repmat(mean(X),N,1))*v1 - (Y-repmat(mean(Y),N,1))*v2;

else % nvar1>nvar2: solve big, joint eigenproblem
    
Sleft  = [zeros(nvar1,nvar1) s12; s21 zeros(nvar2,nvar2)];
Sright = [s11 zeros(nvar1,nvar2); zeros(nvar2,nvar1) s22];
[v,d] = eigen2(Sleft,Sright);
v = fliplr(v);
v1 = v(1:nvar1,1:nvar1);
v2 = v(nvar1+1:end,1:nvar2);   
rho = fliplr(diag(d)'); % highest canonical correlation first
rho(:,nvar2+1:end) = 0;
rho = rho(:,1:nvar1);
    
% scale v1 to give CVs with unit variance
aux1 = v1'*s11*v1; % dispersion of CVs
aux2 = 1./sqrt(diag(aux1));
aux3 = repmat(aux2',nvar1,1);
v1 = v1.*aux3; % now dispersion is unit matrix
%v1'*s11*v1

% sum of correlations between X and CV(X) positive
invstderr = diag(1./std(X));
%invstderrcv = diag(1./sqrt(diag(v1'*s11*v1)));
%sgn = diag(sign(sum(invstderr*s11*v1*invstderrcv)));
sgn = diag(sign(sum(invstderr*s11*v1)));
v1 = v1*sgn;
%figure; bar(sum(invstderr*s11*v1))

% scale v2 to give CVs with unit variance
aux1 = v2'*s22*v2; % dispersion of CVs
aux2 = 1./sqrt(diag(aux1));
aux3 = repmat(aux2',nvar2,1);
v2 = v2.*aux3; % now dispersion is unit matrix
%v2'*s22*v2

% correlations between CV(X) and CV(Y) positive
sgn = diag(sign(diag(v1'*s12*v2)'));
v2 = v2*sgn;
%invstderr = diag(1./std(Y)); figure; bar(sum(invstderr*s22*v2))
%v1'*s12*v2

mads = (X-repmat(mean(X),N,1))*v1 - ...
    [(Y-repmat(mean(Y),N,1))*v2 zeros(N,nvar1-nvar2)];
mads = fliplr(mads);
mads(:,1:nvar1-nvar2) = sqrt(2)*mads(:,1:nvar1-nvar2); % variance = 2
rho = fliplr(rho);

end

%cov(mads)

if nargout>7
    chi2 = sum((mads./repmat(sqrt(2*(1-rho)),N,1)).^2,2); % should be no-change std only
    %chi2 = sum((mads./repmat(std(mads),N,1)).^2,2); % should be no-change std only
    %prob = tanh(1./chi2); % sensible measure but ad hoc
    prob = 1-gammainc(0.5*chi2,0.5*nvar1); % probability of finding larger chi2
end

% output array mads consists of transposed images
mads = reshape(mads,nrows,ncols,nvar1);

return

if nargin==3
    fid = fopen(fnameo,'w');
    fwrite(fid,mads,'float32');
    fclose(fid);
    fid = fopen(strcat(fnameo,'.hdr'),'w'); % write primitive header file
    fprintf(fid,'samples = %d\n',ncols);
    fprintf(fid,'lines   = %d\n',nrows);
    fprintf(fid,'bands   = %d\n',nvar1);
    fprintf(fid,'data type = 4\n');
    fprintf(fid,'rho = {\n'); fprintf(fid,' %g',rho); fprintf(fid,'}\n');
    fclose(fid);
end
