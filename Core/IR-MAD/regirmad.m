function [regirmads,rho,Lambda,v1,v2,s11,s22,s12,w,chi2,CV1,CV2] = regirmad(X,Y,varargin)

% [regirmads,rho,Lambda,v1,v2,s11,s22,s12,w,chi2] = regirmad(X,Y,Lambda,varargin)
%
% Regularized IR-MAD -
%       regularized iteratively reweighted multivariate alteration detection
%
% Input
% X      - file name or 3-D matrix with multivariate band sequential input image
%           number one
% Y      - file name or 3-D matrix with multivariate band sequential input image
%           number two
% Lambda - the regularization or penalization factor,
%          same for both sets of variables (optional, defaults to 0)
% eps    - epsilon for iterations (optional, defaults to 10e-2)
% w      - initial weights (optional)
% 
% Output
% the MAD variates (mads),
% the canonical correlations (rho),
% the regularization or penalization factor, same for both sets of variables (Lambda)
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

% This version gets Lambda for regularization as an input (and outputs it
% again)

% (c) Copyright 2005-2010
% Allan Aasbjerg Nielsen
% aa@space.dtu.dk, www.imm.dtu.dk/~aa
% 8 Sep 2010

maxiter = 30;           % max number of iterations for IR-MAD
graphics = 1;
Lambda = 0;
if nargin>2
    Lambda = varargin{1};
end
eps=10^(-2);            % stop criterion for cancorrs between iterarions for IR-MAD
if nargin>3
    eps = varargin{2};
end
bins = [-inf -6:0.1:6 inf]; % for histxy
%bins = -3:0.1:3; % for histxy
plotbins = bins;
if bins(1) == -inf
    plotbins(1) = bins(2) - (bins(3) - bins(2));
end
if bins(end) == inf
    plotbins(end) = bins(end-1) + (bins(end-1) - bins(end-2));
end

readfile1 = 0; % input1 is 3-D matrix
readfile2 = 0; % input2 is 3-D matrix
if nargin<2, error('Not enough input arguments.'); end
if nargin>5, error('Too many input arguments.'); end
if nargin==5
    w = varargin{end};
    [nroww,ncolw,nvarw] = size(w);
end
if ischar(X), readfile1 = 1; end % input1 is filename
if ischar(Y), readfile2 = 1; end % input2 is filename

if readfile1==1
    [X,p] = freadenvit(X);
    nrows = p(2);
    ncols = p(1);
    nvar1 = p(3);
    X = reshape(X,nrows*ncols,nvar1);
else
    if ndims(X)~=3, error('input1 must be 3-D'); end
    [nrows,ncols,nvar1] = size(X);
    X = reshape(X,nrows*ncols,nvar1);
end

N = nrows*ncols;

if readfile2==1
    [Y,p] = freadenvit(Y);
    nrow2 = p(2);
    ncol2 = p(1);
    nvar2 = p(3);
    Y = reshape(Y,nrow2*ncol2,nvar2);
else
    if ndims(Y)~=3, error('input2 must be 3-D'); end
    [nrow2,ncol2,nvar2] = size(Y);
    Y = reshape(Y,nrow2*ncol2,nvar2);
end

wave = 1:nvar1;

if nvar2~=nvar1
    warning('this implementation requires the same number of variables at the two time points');
end
if nvar2>nvar1
    error('input with highest number of variables must be first set');
end
if ~(ncols==ncol2 && nrows==nrow2)
    error('data in X and Y do not match');
end
if exist('w','var')
    if ~(ncols==ncolw && nrows==nroww)
        error('data in X, Y and w do not match');
    end
else
    w = ones(N,1);
end

% make Omega for penalisation of size, slope and/or curvature in weights
if Lambda ~= 0
%    O1 = eye(nvar1); %size
%    O2 = eye(nvar2);
%    O1 = makeOmega1(nvar1); %slope
%    O2 = makeOmega1(nvar2);
    O1 = makeOmega2(nvar1); %curvature
    O2 = makeOmega2(nvar2);
% same weight on slope and curvature
%    O1 = 3*makeOmega1(nvar1) + makeOmega2(nvar1);
%    O2 = 3*makeOmega1(nvar2) + makeOmega2(nvar2);
end

rho0 = 100*ones(1, nvar1);
rho = [];

for iter=1:maxiter % -------------------------------------------------------

[covxy, meanw] = covw([X Y], w);
meanwX = meanw(1:nvar1);
meanwY = meanw(nvar1+1:end);
% make corr mat from cov mat
% aux = diag(1./sqrt(diag(covxy)));
% covxy = aux * covxy * aux;

s11 = covxy(1:nvar1,1:nvar1);
s22 = covxy(nvar1+1:end,nvar1+1:end);
s12 = covxy(1:nvar1,nvar1+1:end);
s21 = s12';

% change O1 and O2!!!!!!!!!!!!!
if Lambda ~= 0
    O1 = O1*trace(s11)/trace(O1);
    O2 = O2*trace(s22)/trace(O2);
end

A = [zeros(nvar1, nvar2) s12; s21 zeros(nvar2, nvar1)];
%A = [s11 s12; s21 s22];
if Lambda==0
    B = [s11                      zeros(nvar1, nvar2); ...
            zeros(nvar2, nvar1) s22          ];
else
    B = [(1-Lambda)*s11+Lambda*O1 zeros(nvar1, nvar2); ...
            zeros(nvar2, nvar1) (1-Lambda)*s22+Lambda*O2];
end
[V, D, flag3] = eigs(A, B, nvar1, 'LA');
%if ~(flag3==0), warning('*** Convergence problems in eigs ***'); end
rho = [rho; diag(D)'];
%rho = [rho; diag(D)'-1];
V = sqrt(2) * V; % due to normalization in eigs

v1 = V(1:nvar1,:);
v2 = V(nvar1+1:end,:);

% alternative rescaling
%aux1 = v1'*B(1:nvar1,1:nvar1)*v1; % dispersion of CVs
%aux2 = 1./sqrt(diag(aux1));
%aux3 = repmat(aux2',nvar1,1);
%v1 = v1.*aux3; % now dispersion is unit matrix
%%v1'*B(1:nvar1,1:nvar1)*v1
%aux1 = v2'*B(nvar1+1:end,nvar1+1:end)*v2; % dispersion of CVs
%aux2 = 1./sqrt(diag(aux1));
%aux3 = repmat(aux2',nvar2,1);
%v2 = v2.*aux3; % now dispersion is unit matrix
%%v2'*B(nvar1+1:end,nvar1+1:end)*v2

% sum of correlations between X and CV(X) positive
invstderrX = 1./std(X);
%invstderrcv = diag(1./sqrt(diag(v1'*s11*v1)));   % unit matrix here
%sgn = diag(sign(sum(invstderr*s11*v1*invstderrcv)));
sgn = diag(sign(sum(diag(invstderrX)*s11*v1)));
%sgn = diag(sign(sum(s11*v1)));
v1 = v1*sgn;
%figure; bar(sum(invstderrX*s11*v1));
%figure; bar(sum(s11*v1));

% correlations between CV(X) and CV(Y) positive
sgn = diag(sign(diag(v1'*s12*v2)'));
v2 = v2*sgn;
%invstderrY = 1./std(Y); figure; bar(sum(diag(invstderrY)*s22*v2))
%v1'*s12*v2

CV1 = (X-repmat(meanwX,N,1))*v1;
CV2 = (Y-repmat(meanwY,N,1))*v2;

if iter==1
    histxy = hist2d([CV2(:,1) CV1(:,1)], bins, bins);
    figure(100)
    Plot2dHist(histxy, plotbins, plotbins, '', '', '')
    %Plot2dHist(histxy, bins, bins, '', '', '')
    axis equal tight
    %drawnow
    print -depsc2 'histo2dcv.eps'
end

if nvar2<nvar1
    regirmads = CV1 - [zeros(N,nvar1-nvar2) CV2];
%    regirmads = X*v1 - [zeros(N,nvar1-nvar2) Y*v2];
else
    regirmads = CV1 - CV2;
%    regirmads = X*v1 - Y*v2;
end

varmads = 2*(1-rho(end,:));
if Lambda~=0
    varmads = varmads + ...
        Lambda*(sum(v1.*((s11-O1)*v1))+sum(v2.*((s22-O2)*v2)));
end
chi2 = sum((regirmads.^2./repmat(varmads,N,1)),2); % should be no-change std only

if iter==1, disp('Canonical correlations'); end
disp(num2str(rho,' %0.6g'))

if max(abs(rho(end,:)-rho0))<eps break; end

rho0 = rho(end,:);

w = 1-gammainc(0.5*chi2,0.5*nvar1); % probability of finding larger chi2
% Simple spatial extension
%w = imfilter(reshape(w,ncols,nrows),[0 1 0;1 4 1;0 1 0]/8,'symmetric');
%w = imfilter(reshape(w,ncols,nrows),fspecial('gaussian',5,1),'symmetric');

end % end iterations

if iter==maxiter
    warning('No convergence, max number of iterations performed');
end

histxy = hist2d([CV2(:,1) CV1(:,1)], bins, bins);
figure(101)
Plot2dHist(histxy, plotbins, plotbins, '', '', '')
axis equal tight
print -depsc2 'histo2dircv.eps'

figure;
plot(rho, 'o-');
ylabel(strcat('Canonical correlations, \lambda=', num2str(Lambda)));
print -depsc2 'rho1.eps'
print -dtiff 'rho1.tif'
figure
plot(rho', 'o-')
ylabel(strcat('Canonical correlations, \lambda=', num2str(Lambda)));
print -depsc2 'rho2.eps'
print -dtiff 'rho2.tif'

for ii=1:3 % nvar1
    figure;
    subplot(2,1,1)
    plot(wave, v1(:, ii), 'o-');
    ylabel(strcat('Weights for CV1s ',int2str(ii),', \lambda =',num2str(Lambda)));
    subplot(2,1,2)
    plot(wave, v2(:, ii), 'o-');
    ylabel(strcat('Weights for CV2s ',int2str(ii),', \lambda =',num2str(Lambda)));
%    print -depsc2 strcat('w',int2str(ii),'L.eps')
%    print -dtiff strcat('w',int2str(ii),'L.tif')
end

if ~(flag==0), warning('*** Convergence problems in eigs ***'); end

if N<10000, return; end

if graphics
    % no-change obs only
%     aux = diag(1./sqrt(varmads));
%     cor = [diag(invstderrX)*(s11*v1-s12*v2)*aux
%            diag(1./std(Y)) *(s21*v1-s22*v2)*aux];
%     clear aux
    cor = corrcoef([regirmads X Y]); cor = cor(nvar1+1:end, 1:nvar1);
    for ii=1:nvar1
        figure, bar(cor(:,ii)), set(gca, 'ylim', [-1 1]);
        title(strcat(['Correlations between IR-MAD ', int2str(ii), ...
            ' and original variables, \rho = ', num2str(rho(end,ii))]))
    end
    if nvar1>5
        figure
        for ii=1:6
            subplot(3,2,ii), bar(cor(:,ii)), set(gca, 'ylim', [-1 1]);
        end
        print -depsc2 cor1
    end
    if nvar1>11
        figure
        for ii=7:12
            subplot(3,2,ii-6), bar(cor(:,ii)), set(gca, 'ylim', [-1 1]);
        end
        print -depsc2 cor2
    end
%     figure, plotmatrixAA(regirmads)
%     print -deps2 plotmatrix
%     print -djpeg99 plotmatrix
%     set(gcf,'InvertHardCopy','off','Color','k')
%     print -djpeg99 plotmatrixk
    % show MADs with sensible stretching:
    % stretch over -/+ nstd stddevs of no-change observations
    nstd = 20;
    stdpcs = sqrt(varmads);
    % output array consists of transposed images
    regirmads = reshape(regirmads,nrows,ncols,nvar1);
    for iii=1:nvar1
        figure, imshow(regirmads(:,:,iii),nstd*[-stdpcs(iii) stdpcs(iii)])
        title(strcat(['IR-MAD ', int2str(nvar1+1-iii), ', \rho = ', num2str(rho(end,iii))]))
%         print -depsc2 regirmads1
%         print -djpeg  regirmads1
%         set(gcf,'InvertHardCopy','off','Color','k')
%         print -djpeg  regirmads1k
    end
    if nvar1>2
        r = regirmads(:,:,1)/(2*nstd*stdpcs(1))+0.5; %r(r<0)=0; %r(r>1)=1;
        g = regirmads(:,:,2)/(2*nstd*stdpcs(2))+0.5; %g(g<0)=0; %g(g>1)=1;
        b = regirmads(:,:,3)/(2*nstd*stdpcs(3))+0.5; %b(b<0)=0; %b(b>1)=1;
        figure, imshow(cat(3,r,g,b)), %title('IR-MADs 1, 2 and 3 as RGB')
        print -depsc2 regirmadsrgb123
        print -djpeg  regirmadsrgb123
        set(gcf,'InvertHardCopy','off','Color','k')
        print -djpeg  regirmadsrgb123k
    end
%     if nvar1>5
%         r = regirmads(:,:,4)/(2*nstd*stdpcs(4))+0.5; %r(r<0)=0; %r(r>1)=1;
%         g = regirmads(:,:,5)/(2*nstd*stdpcs(5))+0.5; %g(g<0)=0; %g(g>1)=1;
%         b = regirmads(:,:,6)/(2*nstd*stdpcs(6))+0.5; %b(b<0)=0; %b(b>1)=1;
%         figure, imshow(cat(3,r,g,b)), %title('IR-MADs 4, 5 and 6 as RGB')
%         print -depsc2 regirmadsrgb456
%         print -djpeg  regirmadsrgb456
%         set(gcf,'InvertHardCopy','off','Color','k')
%         print -djpeg  regirmadsrgb456k
%     end
% %     if nvar1>6
% %         r = regirmads(:,:,5)/(2*nstd*stdpcs(5))+0.5; %r(r<0)=0; %r(r>1)=1;
% %         g = regirmads(:,:,6)/(2*nstd*stdpcs(6))+0.5; %g(g<0)=0; %g(g>1)=1;
% %         b = regirmads(:,:,7)/(2*nstd*stdpcs(7))+0.5; %b(b<0)=0; %b(b>1)=1;
% %         figure, imshow(cat(3,r,g,b)), %title('IR-MADs 5, 6 and 7 as RGB')
% %         print -depsc2 regirmadsrgb567
% %         print -djpeg  regirmadsrgb567
% %         set(gcf,'InvertHardCopy','off','Color','k')
% %         print -djpeg  regirmadsrgb567k
% %     end
    % imtool(cat(3,r,g,b))
    % imtool(reshape(regirmads(:,:,1)),[], 'InitialMagnification', 100)
    % figure, iii=1; imshow(regirmads(:,:,iii),nstd*[-stdpcs(iii) stdpcs(iii)])
    % iii=1; imtool(regirmads(:,:,iii),nstd*[-stdpcs(iii) stdpcs(iii)], 'InitialMag', 100)
end

% output array mads consists of transposed images
w = reshape(w,nrows,ncols);
chi2 = reshape(chi2,nrows,ncols);
CV1 = reshape(CV1,nrows,ncols,nvar1);
CV2 = reshape(CV2,nrows,ncols,nvar2);

return

fwritehdr(regirmads,'regirmads');
fwritehdr(w,'regirmadsw');
fwritehdr(chi2,'regirmadschi2');

fid = fopen('regirmads.hdr','a');
if fid~=-1
    fprintf(fid, 'rho = {\n');
    for ii = 1:size(rho, 1)
        fprintf(fid,' %g', rho(ii, :)); fprintf(fid, '\n');
    end
    fprintf(fid, '}\n');
end
fclose(fid);
