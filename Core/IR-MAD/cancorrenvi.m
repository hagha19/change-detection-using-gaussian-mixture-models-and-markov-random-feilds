function [CV1,CV2,rho,v1,v2,s11,s22,s12] = cancorrenvi(X,Y)

% [CV1,CV2,rho,v1,v2,s11,s22,s12] = cancorr(X,Y);
%
% CANCORR - canonical correlation analysis
%
% Input
% X      - file name of multivariate band sequential input image
%          number one or image already read into Matlab
% Y      - file name of multivariate band sequential input image
%          number two or image already read into Matlab
% 
% Output
% the canonical variates (CV1 and CV2),
% the canonical correlations (rho),
% the eigenvectors (v1 and v2) normed to give CVs unit variance, and
% the relevant (variance-)covariance matrices.

% (c) Copyright 2005-2010
% Allan Aasbjerg Nielsen, Ph.D., M.Sc.
% aa@space.dtu.dk, www.imm.dtu.dk/~aa
% 19 Sep 2010

if nargin<2, error('Not enough input arguments.'); end
if nargin>2, error('Too many input arguments.'); end
% if nargin==4
%     fnameo1 = varargin{1};
%     if ~ischar(fnameo1), error('fnameo1 should be a char string'); end
%     fnameo2 = varargin{2};
%     if ~ischar(fnameo2), error('fnameo2 should be a char string'); end
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
rho = fliplr(diag(sqrt(d1))'); % highest canonical correlation first
v1 = fliplr(v1);

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
%figure; plot(sum(invstderr*s11*v1),'o')

%%[v2,d2] = eigen2(s21*(s11^(-1))*s12,s22);
%[v2,d2] = eigen2(s21*inv(s11)*s12,s22);
%v2 = fliplr(v2);
%v2 = v2*diag(sign(diag(v2'*s12*v2)));
v2 = invs22*s21*v1; %./repmat(rho,nvar1,1); % scaling doesn't matter
%if nvar2<nvar1
%    v2 = v2(:,1:nvar2);
%    rho(:,end-(nvar1-nvar2)+1:end) = 0;
%end

% scale v2 to give CVs with unit variance
aux1 = v2'*s22*v2; % dispersion of CVs
aux2 = 1./sqrt(diag(aux1));
aux3 = repmat(aux2',nvar2,1);
v2 = v2.*aux3; % now dispersion is unit matrix
%v2'*s22*v2
%invstderr = diag(1./std(Y)); figure; plot(sum(invstderr*s22*v2),'o')

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
%figure; plot(sum(invstderr*s11*v1),'o')

% scale v2 to give CVs with unit variance
aux1 = v2'*s22*v2; % dispersion of CVs
aux2 = 1./sqrt(diag(aux1));
aux3 = repmat(aux2',nvar2,1);
v2 = v2.*aux3; % now dispersion is unit matrix
%v2'*s22*v2

% correlations between CV(X) and CV(Y) positive
sgn = diag(sign(diag(v1'*s12*v2)'));
v2 = v2*sgn;
%invstderr = diag(1./std(Y)); figure; plot(sum(invstderr*s22*v2),'o')
%v1'*s12*v2
    
end

CV1 = (X-repmat(mean(X),N,1))*v1;
CV2 = (Y-repmat(mean(Y),N,1))*v2;
%CV1 = CV1./repmat(std(CV1),N,1); % unit variance
%CV2 = CV2./repmat(std(CV2),N,1); % unit variance

%cov([CV1 CV2])
%MAD = fliplr(CV1-CV2);

% output arrays CV1 and CV2 consist of transposed images
CV1 = reshape(CV1,nrows,ncols,nvar1);
CV2 = reshape(CV2,nrows,ncols,nvar2);

return

if nargin>2
    fid = fopen(fnameo1,'w');
    fwrite(fid,CV1,'float32');
    fclose(fid);
    fid = fopen(fnameo2,'w');
    fwrite(fid,CV2,'float32');
    fclose(fid);
    fid = fopen(strcat(fnameo1,'.hdr'),'w'); % write primitive header file
    fprintf(fid,'samples = %d\n',ncols);
    fprintf(fid,'lines   = %d\n',nrows);
    fprintf(fid,'bands   = %d\n',nvar1);
    fprintf(fid,'data type = 4\n');
    fclose(fid);
    fid = fopen(strcat(fnameo2,'.hdr'),'w'); % write primitive header file
    fprintf(fid,'samples = %d\n',ncols);
    fprintf(fid,'lines   = %d\n',nrows);
    fprintf(fid,'bands   = %d\n',nvar2);
    fprintf(fid,'data type = 4\n');
    fclose(fid);
end
