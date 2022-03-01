function y = spectralEnergy(image, mu, sigma, P, num_cluster)
% obj = gmdistribution(mu,sigma,P);
% Ef = posterior(obj,image);
image=image(:);
mu=mu(:);
sigma=sigma(:);
P=P(:);
for i=1:size(mu,1)
   d = image-mu(i);
   amp = P(i)/sqrt(2*pi*sigma(i));
   y(:,i) = amp*exp(-0.5 * (d.*d)/sigma(i));
end
y = [y(:,1)./(y(:,1)+y(:,2)), y(:,2)./(y(:,1)+y(:,2))];
end