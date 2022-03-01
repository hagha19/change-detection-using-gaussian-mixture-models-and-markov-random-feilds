function comp = ICM1(image, num_cluster, potential, maxItr)
[width,height,bands]=size(image);
image=imstack2vectors(image);
GMmodel = fitgmdist(image, num_cluster);
im_clust = cluster(GMmodel, image);
iter=0;
while(iter<maxItr)
    %     [mu,sigma, P]=GMM_param(image,im_clust,num_cluster);
    mu = GMmodel.mu;
    sigma = GMmodel.Sigma;
    P = GMmodel.PComponents;
    Ef=spectralEnergy(image,mu,sigma, P, num_cluster);
    E1=spatialEnergy(im_clust, width,height,num_cluster);
    E=Ef+potential*E1;
    [tm,segmentation]=min(E,[],2);
    iter=iter+1;
end
comp = zeros(size(im_clust));
b = mean(image(segmentation==1));
a = mean(image(segmentation==2));
[aa, bb] = min([a b]);
if bb==1
comp(segmentation==1)= 2;
comp(segmentation==2)= 1;
else
 comp(segmentation==2)= 2;
comp(segmentation==1)= 1;   
end
comp=reshape(comp,[width height]);

end