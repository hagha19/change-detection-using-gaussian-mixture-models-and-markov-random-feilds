function [im, num_pts] = GMM_RadiometricNormalization(im1 , im2)

[row col band] = size(im1);
dif_im = (single(im1) - single(im2));
reshaped_im = reshape(dif_im , size(im1, 1)* size(im1,2) , size(im1,3));
gmmodel = fitgmdist(reshaped_im, 2);
post = posterior(gmmodel, reshaped_im);
mu = gmmodel.mu;
[~, b] = min(abs(mu));
[id ,~] = find(post(:,mode(b))> .9);
num_pts = length(id);

for ch = 1:band
    a = im1(:,:,ch);
    b = im2(:,:,ch);
    x = b(id);
    t = a(id);
    if 1
    
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
        b_new = feval(fitresult,single(b(:)'));
        im(:,:,ch) = reshape(b_new' , row , col);
        
    end
end
end

    


    
