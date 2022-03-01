function [irmads,rho,v1,v2,s11,s22,s12,w,chi2] = irmadenvi(X,Y,varargin)

% [irmads,rho,v1,v2,s11,s22,s12,w,chi2] = irmadenvi(X,Y,varargin)
%
% IR MAD - iteratively reweighted multivariate alteration detection
%
% Input
% X      - file name of multivariate band sequential input image
%          number one or image already read into Matlab
% Y      - file name of multivariate band sequential input image
%          number two or image already read into Matlab
% epsln  - epsilon for iterations (optional, defaults to 10e-2)
% w      - weights for stats calculation, image (optional)
%
% Output
% the MAD variates (mads),
% the canonical correlations (rho),
% the eigenvectors (v1 and v2) normed to give CVs unit variance,
% the relevant (variance-)covariance matrices, and
% if desired, the final weights applied in calculating the statistics and
% the approximate chi-squared quantity. 
%
% Output array mads consists of transposed images must be viewed with e.g.
%       imshow(reshape(mads(:,:,1),ncols,nrows)',[-3 3])
%
% If disk output is requested a primitive .hdr file for the output file is written;
% the full ENVI .hdr file must be constructed manually.

% To run without iterations use epsln>100

% (c) Copyright 2005-2010
% Allan Aasbjerg Nielsen, Ph.D., M.Sc.
% aa@space.dtu.dk, www.imm.dtu.dk/~aa
% 19 Sep 2010

if nargin<2, error('Not enough input arguments.'); end
if nargin>4, error('Too many input arguments.'); end
maxiter = 30;

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

% if nvar2~=nvar1
%     warning('this implementation requires the same number of variables at the two time points');
% end
if nvar2>nvar1
    error('input with highest number of variables must be first set');
end
if ~(ncols==ncol2 && nrows==nrow2)
    error('data in X and Y do not match');
end

epsln=10^(-2);
if nargin>2
    epsln = varargin{1};
end

if nargin>3
    w = varargin{end};
    [nroww,ncolw,nvarw] = size(w);
end

N = nrows*ncols;
if exist('w','var')
    if ~(ncols==ncolw && nrows==nroww)
        error('data in X, Y and w do not match');
    end
else
    w = ones(N,1);
end

rho0=100*ones(1,nvar1);
rho = [];

for iter=1:maxiter % -------------------------------------------------------
    
[covxy,meanw] = covw([X Y],w);
meanwX = meanw(1:nvar1);
meanwY = meanw(nvar1+1:end);
s11 = covxy(1:nvar1,1:nvar1);
s22 = covxy(nvar1+1:end,nvar1+1:end);
s12 = covxy(1:nvar1,nvar1+1:end);
s21 = s12';

%[v1,d1] = eigen2(s12*(s22^(-1))*s21,s11);
invs22 = inv(s22);
[v1,d1] = eigen2(s12*invs22*s21,s11);
rho = [rho; diag(sqrt(d1))']; % lowest canonical correlation first

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

%%[v2,d2] = eigen2(s21*(s11^(-1))*s12,s22);
%[v2,d2] = eigen2(s21*inv(s11)*s12,s22);
%v2 = v2*diag(sign(diag(v1'*s12*v2)));
v2 = invs22*s21*v1; %./repmat(rho(end,:),nvar1,1); % scaling doesn't matter

% scale v2 to give CVs with unit variance
aux1 = v2'*s22*v2; % dispersion of CVs
aux2 = 1./sqrt(diag(aux1));
aux3 = repmat(aux2',nvar2,1);
v2 = v2.*aux3; % now dispersion is unit matrix
%v2'*s22*v2
if nvar2<nvar1
    v2 = v2(:,nvar1-nvar2+1:end);
    rho(:,1:nvar1-nvar2) = 0;
end

if nvar2<nvar1
    irmads = (X-repmat(meanwX,N,1))*v1 - ...
        [zeros(N,nvar1-nvar2) (Y-repmat(meanwY,N,1))*v2];
    % variance = 2 = 2*(1-rho), rho=0
    irmads(:,1:nvar1-nvar2) = sqrt(2)*irmads(:,1:nvar1-nvar2);
else
    irmads = (X-repmat(meanwX,N,1))*v1 - (Y-repmat(meanwY,N,1))*v2;
end

if iter==1, disp('Canonical correlations'); end
disp(num2str(rho(end,:),' %0.6g'))

varmads = 2*(1-rho(end,:));
chi2 = sum((irmads.^2./repmat(varmads,N,1)),2); % should be no-change std only

if max(abs(rho(end,:)-rho0))<epsln break; end
rho0 = rho(end,:);
w = 1-gammainc(0.5*chi2,0.5*nvar1); % probability of finding larger chi2

end % end iterations

if iter==maxiter
    warning('No convergence, max number of iterations performed');
end

irmads = reshape(irmads,nrows,ncols,nvar1);

if nargout>7
    w = reshape(w,nrows,ncols);
end
if nargout>8
    chi2 = reshape(chi2,nrows,ncols);
end

while 0

% show MADs with sensible stretching:
% stretch over -/+ nstd stddevs of no-change observations
nstd = 10;
stdmad = sqrt(2*(1-rho(end,:)));

% MAD corresponding to highest order CVs
figure, imshow(irmads(:,:,end),nstd*[-stdmad(end) stdmad(end)])
% MAD corresponding to lowest order CVs
figure, imshow(irmads(:,:,1),nstd*[-stdmad(1) stdmad(1)])

% MADs corresponding to highest order CVs, RGB
r = irmads(:,:,end  )/(2*nstd*stdmad(end  ))+0.5; r(r<0)=0; r(r>1)=1;
g = irmads(:,:,end-1)/(2*nstd*stdmad(end-1))+0.5; g(g<0)=0; g(g>1)=1;
b = irmads(:,:,end-2)/(2*nstd*stdmad(end-2))+0.5; b(b<0)=0; b(b>1)=1;
figure, imshow(cat(3,r, g, b))
%colormap(rgb);

% MADs corresponding to lowest order CVs, RGB
r = irmads(:,:,1)/(2*nstd*stdmad(1))+0.5; r(r<0)=0; r(r>1)=1;
g = irmads(:,:,2)/(2*nstd*stdmad(2))+0.5; g(g<0)=0; g(g>1)=1;
b = irmads(:,:,3)/(2*nstd*stdmad(3))+0.5; b(b<0)=0; b(b>1)=1;
figure, imshow(cat(3,r, g, b))
%colormap(rgb);

end % while 0

return

% output like this or use 'fwritehdr'
if nargin==4
    fid = fopen('irmads','w');
    fwrite(fid,irmads,'float32');
    fclose(fid);
    fid = fopen(strcat('irmads','.hdr'),'w'); % write primitive header file
    fprintf(fid,'samples = %d\n',ncols);
    fprintf(fid,'lines   = %d\n',nrows);
    fprintf(fid,'bands   = %d\n',nvar1);
    fprintf(fid,'data type = 4\n');
    fprintf(fid,'epsilon = %g\n',epsln);
    fprintf(fid,'iterations = %d\n',iter);
%     fprintf(fid,'rho = {\n'); fprintf(fid,' %g',rho); fprintf(fid,'}\n');
    fid = fopen('irmads.hdr','a');
    if fid~=-1
        fprintf(fid, 'rho = {\n');
        for ii = 1:size(rho, 1)
            fprintf(fid,' %g', rho(ii, :)); fprintf(fid, '\n');
        end
        fprintf(fid, '}\n');
    end
    fclose(fid);
end
