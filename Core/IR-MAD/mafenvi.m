function [mafs,ac,v,d,sigmad,sigma] = mafenvi(X,varargin)

% [mafs,ac,v,d,sigmad,sigma] = mafenvi(X,varargin)
%
% MAF - maximum autocorrelation factor analysis
%
% Input
% X      - file name of multivariate band sequential input image or
%          variable name
% flag   - 0, 1 or 2 (optional, defaults to 1)
%
% flag = 0 - eigenvectors scaled as from eig
% flag = 1 - MAFs will be standardised to unit variance (default)
% flag = 2 - eigenvectors will be unit vectors
% flag = 3 - noise in model has unit variance
%
% Output
% mafs   - the factors stored in an ncols by nrows by nvars array
%          (each factor is transposed, see below)
% ac     - autocorrelation in each factor
% v      - the eigenvectors or weights to obtain the factors from the
%          original variables
% d      - the eigenvalues
% sigmad - variance-covariance matrix of difference between spatially
%          shifted and original images
% sigma  - variance-covariance matrix of input image

% (c) Copyright 2005-2010
% Allan Aasbjerg Nielsen, Ph.D., M.Sc.
% aa@space.dtu.dk, www.imm.dtu.dk/~aa
% 19 Sep 2010

if nargin<1, error('Not enough input arguments.'); end
if nargin>2, error('Too many input arguments.'); end
readfile1 = 0; % input1 is 3-D matrix

if ischar(X), readfile1 = 1; end % input1 is filename

flag = 1;
if nargin>1
    flag = varargin{1};
end

if readfile1==1
    [X,p] = freadenvit(X);
    nrows = p(2); ncols = p(1); nvars = p(3);
else
    if ndims(X)~=3, error('X must be 3-D'); end
    [nrows,ncols,nvars] = size(X);
end
X = reshape(X,nrows*ncols,nvars);

sigma = cov(X);
sigmad = pool(X,nrows,ncols,nvars);

[v,d] = eigen2(sigmad,sigma);
d=diag(d)';
ac = 1-0.5*d; % autocorrelation

% sum of correlations between X and MAFs positive
invstderr = diag(1./std(X));
invstderrmaf = diag(1./sqrt(diag(v'*sigma*v)));
sgn = diag(sign(sum(invstderr*sigma*v*invstderrmaf)));
v = v*sgn;

N = nrows*ncols;
X = X-repmat(mean(X),N,1);
if flag==1 % unit variance
    % scale v to give MAFs with unit variance
    aux1 = v'*sigma*v; % dispersion of MAFs
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
elseif flag==3 % noise part has unit variance (SMAF)
    aux1 = 0.5*v'*sigmad*v; % sigma_n = sigma_d/2
    aux2 = 1./sqrt(diag(aux1));
    aux3 = repmat(aux2',nvars,1);
    v = v.*aux3; % now noise dispersion is unit matrix
    %0.5*v'*sigmad*v
end
mafs = X*v;
%cov(mafs)

% output array mafs consists of transposed images
mafs = reshape(mafs,nrows,ncols,nvars);

return

if nargin==3
    fid = fopen(fnameo,'w');
    fwrite(fid,mafs,'float32');
    fclose(fid);
    fid = fopen(strcat(fnameo,'.hdr'),'w'); % write primitive header file
    fprintf(fid,'samples = %d\n',ncols);
    fprintf(fid,'lines   = %d\n',nrows);
    fprintf(fid,'bands   = %d\n',nvars);
    fprintf(fid,'data type = 4\n');
    fclose(fid);
end
