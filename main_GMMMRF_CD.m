function varargout = main_GMMMRF_CD(varargin)
% MAIN_GMMMRF_CD MATLAB code for main_GMMMRF_CD.fig
%      MAIN_GMMMRF_CD, by itself, creates a new MAIN_GMMMRF_CD or raises the existing
%      singleton*.
%
%      H = MAIN_GMMMRF_CD returns the handle to a new MAIN_GMMMRF_CD or the handle to
%      the existing singleton*.
%
%      MAIN_GMMMRF_CD('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MAIN_GMMMRF_CD.M with the given input arguments.
%
%      MAIN_GMMMRF_CD('Property','Value',...) creates a new MAIN_GMMMRF_CD or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before main_GMMMRF_CD_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to main_GMMMRF_CD_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help main_GMMMRF_CD

% Last Modified by GUIDE v2.5 04-Nov-2017 20:37:19

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @main_GMMMRF_CD_OpeningFcn, ...
                   'gui_OutputFcn',  @main_GMMMRF_CD_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before main_GMMMRF_CD is made visible.
function main_GMMMRF_CD_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to main_GMMMRF_CD (see VARARGIN)

% Choose default command line output for main_GMMMRF_CD
handles.output = hObject;
addpath(genpath('Core'));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes main_GMMMRF_CD wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = main_GMMMRF_CD_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
warning off
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
a = subGUI();
pathname1 = a.pathname1;
handles.pathname1 = pathname1;
pathname2 = a.pathname2;
handles.pathname2 = pathname2;
handles.TransType = a.TransType;
handles.RadiometricMethod = a.RadiometricMethod;
handles.saveChangeMap = a.saveChangeMap;
handles.saveNormsecondIm = a.saveNormsecondIm;
handles.saveLog = a.saveLog;
id = strfind(pathname1, '\');
txtpath = [pathname1(1:id(end)),'Log.txt'];
if handles.saveLog
    dlmwrite(txtpath, 'Automatic Change Detection Analysis of the Difference Image using Gaussian mixture model and Markov random fields',...
        'delimiter','');
    dlmwrite(txtpath , '%-----University of Tehran','-append','delimiter','');
    dlmwrite(txtpath , '%----Programmed By Hamid Ghanbari','-append','delimiter','');
    dlmwrite(txtpath , ['first image path: ',pathname1],'-append','delimiter','');
    dlmwrite(txtpath , ['second image path: ',pathname2],'-append','delimiter','');
    if handles.TransType==1
    dlmwrite(txtpath , ['RCSS Transformation Type: ','First Order PolyNomial'],'-append','delimiter','');
    else
      dlmwrite(txtpath , ['RCSS Transformation Type: ','Artificial Neural Network'],'-append','delimiter','');
    end
    if handles.RadiometricMethod==1
    dlmwrite(txtpath , ['RAdiometric Normalization Approach: ','IR-MAD'],'-append','delimiter','');
    else
      dlmwrite(txtpath , ['RAdiometric Normalization Approach: ','GMM based'],'-append','delimiter','');
    end
    
end
handles.txtpath = txtpath;
guidata(hObject, handles)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
saveNormsecondIm = handles.saveNormsecondIm;
saveLog = handles.saveLog;
txtpath = handles.txtpath;
pathname1 = handles.pathname1;
pathname2 = handles.pathname2;
TransType = handles.TransType;
RadiometricMethod = handles.RadiometricMethod;
im1 = imread(pathname1);
im2 = imread(pathname2);
[m,n,o] = size(im1);
[m1,n1,o1] = size(im2);
if m1~= m | n1~=n | o1~=o
    warndlg('the size of input images is not the same: please select another images','!! image dimension warning !!')
    Return;
end
if saveLog
    dlmwrite(txtpath , ['Image Length: ',num2str(m)],'-append','delimiter','');
    dlmwrite(txtpath , ['Image Width: ',num2str(n)],'-append','delimiter','');
    dlmwrite(txtpath , ['SamplesPerPixel: ',num2str(o)],'-append','delimiter','');
    
end
classType = {'uint16' , 'uint8'};
class1 = class(im1);
class2 = class(im2);
if RadiometricMethod==1 & size(im1,3)==3
[im2, num_pts] = RadiometricNormalization(single(im1) , single(im2));
else
[im2, num_pts] = GMM_RadiometricNormalization(single(im1) , single(im2));
end

if saveLog
    dlmwrite(txtpath , ['number of radiometric control set samples : ',num2str(num_pts)],'-append','delimiter','');
    
end


if find(strcmpi(class1,classType ) , 1)-1
    im1 = uint8(im1);
else
    im1 = uint16(im1);
end
if find(strcmpi(class2,classType ) , 1)-1
    im2 = uint8(im2);
    if saveNormsecondIm
        id = strfind(txtpath, '\');
        filename = [txtpath(1:id(end)),'NormalizedSecondImage.tif'];
        try
        multibandwritetiff(im2,filename,8);
        end
    end
else
    im2 = uint16(im2);
    if saveNormsecondIm
        id = strfind(txtpath, '\');
        filename = [txtpath(1:id(end)),'NormalizedSecondImage.tif'];
        try
        multibandwritetiff(im2,filename,16);
        end
    end
end
set(handles.axes1, 'visible', 'on');
set(handles.axes2, 'visible', 'on');
set(handles.axes1, 'xTick',[], 'yTick', []);
set(handles.axes2, 'xTick',[], 'yTick', []);

if size(im1,3)==1
     axes(handles.axes1);
imshow(im1(:,:,1),[]);
axes(handles.axes2);
imshow(im2(:,:,1),[]);

else
axes(handles.axes1);
imshow(im1(:,:,1:3),[]);
axes(handles.axes2);
imshow(im2(:,:,1:3),[]);
end
set(handles.text1, 'visible', 'on');
 set(handles.text1, 'string', 'First Image');
 set(handles.text2, 'visible', 'on');
set(handles.text2, 'string', 'Second Image');
handles.im1 = im1;
handles.im2 = im2;
handles.class1 = class1;
handles.class2 = class2;
guidata(hObject, handles);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
saveChangeMap = handles.saveChangeMap;
txtpath = handles.txtpath;
im1 = handles.im1;
im2 = handles.im2;
im1 = single(im1);
im2 = single(im2);
%% EXtracting image diferences
img_dif = abs(im1 - im2);
% img_ratioing =  double(im1)./double(im2);

imge_difference = sqrt(sum(img_dif.^2 , 3));
% imge_difference = sqrt(sum(abs(img_dif) , 3));

max_diff = max(imge_difference(:));
min_diff = min(imge_difference(:));
imge_difference = round(((imge_difference-min_diff)/max_diff)*65535);
% imge_difference(:,:,6:10) = img_ratioing;
class_num =2;
potential = 0.5;
maxItr = 20;
%% diference image clustering
segmentation = ICM1(imge_difference,class_num,potential,maxItr);
handles.segmentation= segmentation;
% global im1 im2 segmentation
class1 = handles.class1;
class2 = handles.class2;
classType = {'uint16' , 'uint8'};
if find(strcmpi(class1,classType ) , 1) - 1
    im1 = uint8(im1);
else
    im1 = uint16(im1);
end
if find(strcmpi(class2,classType ) , 1)-1
    im2 = uint8(im2);
else
    im2 = uint16(im2);
end
% save change map
if saveChangeMap
segmentation = segmentation-1;
id = strfind(txtpath, '\');
filename = [txtpath(1:id(end)),'ChangeMap.tif'];
imwrite(im2uint8(segmentation),filename);    
end

set(handles.axes3, 'visible', 'on');
set(handles.axes3, 'xTick', [],'yTick', []);
set(handles.axes4, 'visible', 'on');
set(handles.axes4, 'xTick', [],'yTick', []);
axes(handles.axes3);
imshow(segmentation,[])
zoomerf(handles, im1, im2, segmentation);
guidata(hObject, handles)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
pathname1 = handles.pathname1;
segmentation = handles.segmentation;
segmentation = segmentation-1;
[pathstr,name,ext] = fileparts(pathname1);
imwrite(im2uint8(segmentation),[pathstr,'\ChangeMap.tif']);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pushbutton2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
