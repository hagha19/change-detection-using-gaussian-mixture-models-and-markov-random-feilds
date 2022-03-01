function img = zoomerf(handles, im1, im2, segmentation)
axes(handles.axes3)
imshow(segmentation,[]);
initialPosition=[10 10 80 80];
im1Crop = imrect(handles.axes3, initialPosition);
im1Crop.addNewPositionCallback (@(pos)mycallback(pos, handles, im1, im2, segmentation));
end


function mycallback(pos, handles, im1, im2, segmentation)
x1 = pos(1);
x2 = pos(1)+pos(3);
y1 = pos(2);
y2 = pos(2)+pos(4);
im11 = im1(round(y1:y2) , round(x1:x2) , :);
im22 = im2(round(y1:y2) , round(x1:x2) , :);
seg = segmentation(round(y1:y2) , round(x1:x2) , :);
axes(handles.axes1)
imshow(im11(:,:,1:3), []);
axes(handles.axes2)
imshow(im22(:,:,1:3), []);
axes(handles.axes4)
imshow(seg,[]);

global pos
end