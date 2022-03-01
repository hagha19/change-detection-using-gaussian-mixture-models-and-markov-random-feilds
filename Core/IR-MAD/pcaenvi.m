function [pcs,v,d,sigma] = pcaenvi(X,varargin)

% [pcs,v,d,sigma] = pcaenvi(X,varargin)
%
% PCA - principal component analysis
%
% Input
% X      - file name of multivariate band sequential input image or
%          variable name
% flag   - 0, 1 or 2 (optional, defaults to 2)
%
% flag = 0 - eigenvectors scaled as from eig
% flag = 1 - PCs will be standardised to unit variance
% flag = 2 - eigenvectors will be unit vectors (default)
%
% Output
% pcs    - the factors stored in an ncols by nrows by nvars array
%          (each factor is transposed, see below)
% v      - the eigenvectors or weights to obtain the factors from the
%          original variables
% d      - the eigenvalues
% sigma  - variance-covariance matrix of input image

% (c) Copyright 2005-2010
% Allan Aasbjerg Nielsen, Ph.D., M.Sc.
% aa@space.dtu.dk, www.imm.dtu.dk/~aa
% 19 Jan 2010

if nargin<1, error('Not enough input arguments.'); end
if nargin>2, error('Too many input arguments.'); end
readfile1 = 0; % input1 is 3-D matrix

if ischar(X), readfile1 = 1; end % input1 is filename

flag = 2;
if nargin>1
    flag = varargin{1};
end
% if nargin==3
%     fnameo = varargin{end};
%     if ~ischar(fnameo), error('fnameo should be a char string'); end
% end

if readfile1==1
    [X,p] = freadenvit(X);
    nrows = p(2); ncols = p(1); nvars = p(3);
else
    if ndims(X)~=3, error('X must be 3-D'); end
    [nrows, ncols, nvars] = size(X);
end
X = reshape(X,nrows*ncols,nvars);

sigma = cov(X);

[v1,d1] = eig(sigma);
d2 = diag(d1);
[dum,idx] = sort(d2);
d = fliplr(d2(idx)'); % largest eigenvalue first
v = fliplr(v1(:,idx));

% sum of correlations between X and PCs positive
invstderr = diag(1./std(X));
invstderrpc = diag(1./sqrt(diag(v'*sigma*v)));
sgn = diag(sign(sum(invstderr*sigma*v*invstderrpc)));
v = v*sgn;

N = nrows*ncols;
X = X-repmat(mean(X),N,1);
if flag==1 % unit variance
    % scale v to give PCs with unit variance
    aux1 = v'*sigma*v; % dispersion of PCs
    aux2 = 1./sqrt(diag(aux1));
    aux3 = repmat(aux2',nvars,1);
    v = v.*aux3; % now dispersion is unit matrix
    %v'*sigma*v
elseif flag==2 % v unit vector
    aux1 = v'*v;
    aux2 = 1./sqrt(diag(aux1));
    aux3 = repmat(aux2',nvars,1);
    v = v.*aux3; % now v contains unit vectors in columns
    %v'*v
end
pcs = X*v;
%cov(pcs)

% output array pcs consists of transposed images
pcs = reshape(pcs,nrows,ncols,nvars);

return

if nargin==3
    fid = fopen(fnameo,'w');
    fwrite(fid,pcs,'float32');
    fclose(fid);
    fid = fopen(strcat(fnameo,'.hdr'),'w'); % write primitive header file
    fprintf(fid,'samples = %d\n',ncols);
    fprintf(fid,'lines   = %d\n',nrows);
    fprintf(fid,'bands   = %d\n',nvars);
    fprintf(fid,'data type = 4\n');
    fclose(fid);
end
