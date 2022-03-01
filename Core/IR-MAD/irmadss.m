% Iteratively re-weighted MAD in scale space

% (c) Copyright 2007-2010
% Allan Aasbjerg Nielsen, Ph.D., M.Sc.
% aa@space.dtu.dk, www.imm.dtu.dk/~aa
% 19 Sep 2010

X0 = freadenvit('C:\cygwin\home\aa\images\20010626bsq');
Y0 = freadenvit('C:\cygwin\home\aa\images\20010829bsq');
% X0 = freadenvit('C:\cygwin\home\aa\ENVI\reservoirt1');
% Y0 = freadenvit('C:\cygwin\home\aa\ENVI\reservoirt2');
% X0 = freadenvit('C:\cygwin\home\aa\brisbaneT1');
% Y0 = freadenvit('C:\cygwin\home\aa\brisbaneT2');
% X0 = freadenvit('C:\cygwin\home\image\images\hips\thika87');
% Y0 = freadenvit('C:\cygwin\home\image\images\hips\thika89');
% X0 = freadenvit('C:\cygwin\home\image\images\hips\avirisBand');
% Y0 = X0(:,:,1:15);
% X0 = X0(:,:,16:30);
% images too big
% X0 = freadenvit('C:\cygwin\home\aa\images\Afghanistan\Spot\20070314_subset');
% Y0 = freadenvit('C:\cygwin\home\aa\images\Afghanistan\Spot\20070505_subset');

%h = [0.05 0.25 0.40 0.25 0.05]'; h2 = h*h';
h2 = fspecial('gaussian',5,1);
[ncols nrows nvars] = size(X0);
w = ones(ncols, nrows);
rhos = [];
nres = 5; % number of scales including full resolution
epsi = 0.01;

% Perform analysis in scale space on smoothed versions of data,
% no sub-sampling
for res=1:nres-1
    X = X0;
    Y = Y0;
    for rres=1:nres-res
        %res, rres
        X = imfilter(X, h2, 'symmetric');
        Y = imfilter(Y, h2, 'symmetric');
    end
    [irmads, rho, v1, v2, s11, s22, s12, w, chi2] = irmadenvi(X, Y, epsi, w);
    figure
    %imshow(X(:,:,4), [])
    imshow(w, [0 1])
    %imshow(chi2, [])
    drawnow
    %truesize
    rhos = [rhos; rho; NaN*ones(1,nvars)];
    clear irmads rho v1 v2 s22 s22 s12 chi2;
    %keyboard
end
[irmads, rho, v1, v2, s11, s22, s12, w, chi2] = irmadenvi(X0, Y0, epsi, w);
figure
%imshow(X(:,:,4), [])
imshow(w, [0 1])
%imshow(chi2',[])
drawnow
%truesize
rhos = [rhos; rho];
figure, plot(rhos, 'o-')
%set(gcf, 'InvertHardCopy', 'on')

while 0

% show MADs with sensible stretching:
% stretch over -/+ nstd stddevs of no-change observations
nstd = 10;
stdmad = sqrt(2*(1-rho(end, :)));

% MAD corresponding to highest order CVs
figure, imshow(irmads(:,:,end)',nstd*[-stdmad(end) stdmad(end)])
%truesize
% MAD corresponding to lowest order CVs
figure, imshow(irmads(:, :, 1),nstd*[-stdmad(1)   stdmad(1)])
%truesize

% MADs corresponding to highest order CVs
r = irmads(:,:,end  )/(2*nstd*stdmad(end  ))+0.5; %r(r<0)=0; r(r>1)=1;
g = irmads(:,:,end-1)/(2*nstd*stdmad(end-1))+0.5; %g(g<0)=0; g(g>1)=1;
b = irmads(:,:,end-2)/(2*nstd*stdmad(end-2))+0.5; %b(b<0)=0; b(b>1)=1;
figure, imshow(cat(3, r, g, b))
%colormap(rgb);
%truesize

% MADs corresponding to lowest order CVs
r = irmads(:,:,1)/(2*nstd*stdmad(1))+0.5; %r(r<0)=0; r(r>1)=1;
g = irmads(:,:,2)/(2*nstd*stdmad(2))+0.5; %g(g<0)=0; g(g>1)=1;
b = irmads(:,:,3)/(2*nstd*stdmad(3))+0.5; %b(b<0)=0; b(b>1)=1;
figure, imshow(cat(3, r, g, b))
%colormap(rgb);
%truesize

end % while 0

while 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform analysis in scale space on smoothed versions of data
% one scale at a time
[ncols nrows nvars] = size(X0);
w = ones(ncols, nrows);
for rres=2:nres
    X = imfilter(X, h2, 'symmetric');
    Y = imfilter(Y, h2, 'symmetric');
end
[irmads, rho, v1, v2, s11, s22, s12, w, chi2] = irmadenvi(X0, Y0, epsi, w);
figure
%imshow(X(:,:,4), [])
imshow(w, [0 1])
%imshow(chi2, [])
drawnow
%truesize
end % while 0
