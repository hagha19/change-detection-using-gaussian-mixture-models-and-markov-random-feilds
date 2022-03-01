function [im, num_pts] = RadiometricNormalization(im1 , im2)
% this function performs radiometric normalization using IR-MAD algorithm 
% INPUTS:
% im1: first image acquired in the first time
% im2: second image acquired in the second time
% OUTPUT:
% im: second image which has changed in order to has equal radiometric
% conditions to first image.
mod = 0;
addpath(genpath('imm4695'));
[mads,rho,v1,v2,s11,s22,s12,prob,chi2] = madenvi(im1,im2);
sum_mad = 0;
[row col band] = size(im1);
for ch = 1:band
    sigma(ch) = sqrt(var(var(mads(:,:,ch))));
    sum_mad = (mads(:,:,ch)/sigma(ch)).^2+sum_mad;
end
T = chi2inv(.99 , band);
id = find(sum_mad<T);
num_pts = length(id);

for ch = 1:band
    a = im1(:,:,ch);
    b = im2(:,:,ch);
    x = b(id);
    t = a(id);
    if mod
    
    % Choose a Training Function
    % For a list of all training functions type: help nntrain
    % 'trainlm' is usually fastest.
    % 'trainbr' takes longer but may be better for challenging problems.
    % 'trainscg' uses less memory. NFTOOL falls back to this in low memory situations.
    trainFcn = 'trainlm';  % Levenberg-Marquardt
    
    % Create a Fitting Network
    hiddenLayerSize = 1;
    net = fitnet(hiddenLayerSize,trainFcn);
    
    % Setup Division of Data for Training, Validation, Testing
    net.divideParam.trainRatio = 70/100;
    net.divideParam.valRatio = 15/100;
    net.divideParam.testRatio = 15/100;
    net.trainParam.epochs = 100;
    net.trainParam.showWindow = 0;
    % Train the Network
    [net,tr] = train(net,x',t');
%     view(net)
    % Test the Network
    b_new = net(b(:)');
    im(:,:,ch) = reshape(b_new' , row , col);
    else
        [xData, yData] = prepareCurveData( x, t );
        
        % Set up fittype and options.
        ft = fittype( 'poly1' );
        
        % Fit model to data.
        [fitresult, ~] = fit( xData, yData, ft );
        b_new = feval(fitresult,b(:)');
        im(:,:,ch) = reshape(b_new' , row , col);
        
    end
end
end

    


    
