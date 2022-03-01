function multibandwritetiff(im,filename,bits)
t = Tiff(filename, 'w')
t.setTag('Photometric',Tiff.Photometric.MinIsBlack );
t.setTag('ImageLength',size(im,1));
t.setTag('ImageWidth',size(im,2));
t.setTag('BitsPerSample',bits);
t.setTag('SamplesPerPixel',size(im,3));
t.setTag('ExtraSamples',Tiff.ExtraSamples.UnassociatedAlpha);
t.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
t.write(im)
end