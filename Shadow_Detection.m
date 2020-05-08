function varargout = Shadow_Detection(varargin)
% SHADOW_DETECTION MATLAB code for Shadow_Detection.fig
%      SHADOW_DETECTION, by itself, creates a new SHADOW_DETECTION or raises the existing
%      singleton*.
%
%      H = SHADOW_DETECTION returns the handle to a new SHADOW_DETECTION or the handle to
%      the existing singleton*.
%
%      SHADOW_DETECTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHADOW_DETECTION.M with the given input arguments.
%
%      SHADOW_DETECTION('Property','Value',...) creates a new SHADOW_DETECTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Shadow_Detection_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Shadow_Detection_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Shadow_Detection

% Last Modified by GUIDE v2.5 14-Apr-2019 14:23:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Shadow_Detection_OpeningFcn, ...
                   'gui_OutputFcn',  @Shadow_Detection_OutputFcn, ...
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


% --- Executes just before Shadow_Detection is made visible.
function Shadow_Detection_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Shadow_Detection (see VARARGIN)

% Choose default command line output for Shadow_Detection
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Shadow_Detection wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Shadow_Detection_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Image 
[filename, pathname] = uigetfile({'*.*'},'File Selector');
 A = strcat(pathname, filename);
 Image = imread(A);
 axes(handles.axes1);
 imshow(A)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Image
I=rgb2gray(Image);
%%=================================================================================================
[r,c]=size(I);
mark=zeros(r,c);

n=zeros(1,256);
B=zeros(r,c);
threshold1=0;
preset=1;
x=preset+1;
while(x>preset)
    
% Compute the histogram============================================================================
for l=1:r
 for m=1:c
     if(mark(l,m)==0)
        z=I(l,m);
        n(z+1)=n(z+1)+1;
     end
end
end
%============================================================================
N=sum(n); % sum the values of all the histogram values
max=0; %initialize maximum to zero
%%================================================================================================
for l=1:256
    P(l)=n(l)/N; %Computing the probability of each intensity level
end
%%================================================================================================
threshold2=threshold1;
for T1=2:255      % step through all thresholds from 2 to 255
    w0=sum(P(1:T1)); % Probability of class 1 (separated by threshold)
    w1=sum(P(T1+1:256)); %probability of class2 (separated by threshold)
    u0=dot([0:T1-1],P(1:T1))/w0; % class mean u0
    u1=dot([T1:255],P(T1+1:256))/w1; % class mean u1
    sigma=w0*w1*((u1-u0)^2); % compute sigma i.e variance(between class)
    if sigma>max % compare sigma with maximum 
        max=sigma; % update the value of max i.e max=sigma
        threshold1=T1-1; % desired threshold corresponds to maximum variance of between class
    end
end
x=abs(threshold1-threshold2);
%%====================================================================================================
H=zeros(r,c);
for l=1:r
 for m=1:c
     if(mark(l,m)==0)
        H(l,m) = I(l,m);
     end
end
end
 %Divide the image into two using threshold
  R0=H (H<threshold1);
  R1=H (H>=threshold1);
  %calculate average of each segment
  miu0=mean(R0(:));
  miu1=mean(R1(:));
 
%%=======================segmenting the image ========================================================
for l=1:r
    for m=1:c
        if(mark(l,m)==0)
            if(I(l,m)<miu0)
                B(l,m)=miu0;
                mark(l,m)=1;
            end
            if(I(l,m)>miu1)
                B(l,m)=miu1;
                mark(l,m)=1;
            end
        end
    end
end
end
% Display the converted Image
axes(handles.axes2);
imshow(B)
% ====================================================================================================



% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Image
x = imadjust(Image,[.2 .3 0; .6 .7 1],[],[]);
%converting to HSV color model
hsv = rgb2hsv(x);
h = medfilt2(hsv(:,:,1),[5 5]);  %median filtering the hue component to remove the noise
s = (hsv(:,:,2));
v = (hsv(:,:,3));
%creating the shadow ratio using color invariant model
shadow_ratio = ((4/pi).*atan(((v))./((s)+((h)))));
%creating shadow mask
%shadow_mask1 = imbinarize(shadow_ratio,'adaptive'); %adaptive thresholding
shadow_mask2 = shadow_ratio>0.5; %global thresholding
shadow_mask3 = shadow_ratio>0.3; 
shadow_mask4 = shadow_ratio>0.7; 
%creating the histogram of the shadow_ratio
axes(handles.axes3);
imshow(shadow_mask2)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Image
b = histeq(rgb2gray(Image),64);  %%%histogram equilization
ima1=double(b);   %convert the image to double
%ima1 = ima1(:,:,3);
%%%%ima11=ima1(:); %converting image matrix to vector
%%%copy=ima1;  %copy the image
[R,C]=size(ima1); %to get the dimension of the matrix
m=R*C;
co1=zeros(m,1);
co2=zeros(m,1);
co3=zeros(m,1);
%make the centroid matrix
mu=[20 117 240]; %initialization of the matrix
%recalculation of the centroid
ass1=zeros(R,C); %zeroes matrix for the first cluster
ass2=zeros(R,C); %zeroes matrix for the second cluster
ass3=zeros(R,C); %zeroes matrix for the third cluster
while(true)
    oldmu=mu;  %store the previous centroid values for comparison later
    for i=1:R %row
        for j=1:C  %column
            r=ima1(i,j); %pixel value of the i,j coordinate
            ab=abs(sqrt((ima1(i,j)^2)-mu.^2)); %to find out that the pixel will  belong to which cluster
            mn=find(ab==min(ab));
            if mn(1)==1
                ass1(i,j)=r;  %assigning to the first cluster
            end
            if mn(1)==2
                ass2(i,j)=r;   %assigning to the second cluster
            end
            if mn(1)==3
                ass3(i,j)=r;   %assigning to the third cluster
            end
        end
    end
   %finding the average of the clusters to determine the new centroid of
   %the clusters
    co1=ass1(:);  %transfer into vector
    su1=sum(co1);  %sum of the vector
    fi1=find(co1);  %to find non-zero element
    len1=length(fi1); %to find the length of non-zero element
    mm1=su1/len1;
    mm11=round(mm1);  %new centre element
    %now to calculate second element of centroid
    co2=ass2(:);
    su2=sum(co2);
    fi2=find(co2);
    len2=length(fi2);
    mm2=su2/len2;
    mm22=round(mm2);
    %now to calculate third elemnt of centroid
    co3=ass3(:);
    su3=sum(co3);
    fi3=find(co3);
    len3=length(fi3);
    mm3=su3/len3;
    mm33=round(mm3);
    %new centroid
    mu=[mm11 mm22 mm33];
    if(isequal(mu, oldmu))  %%%%changed here because of error 
        break;  
    end
 end
%labelling of the clusters
for i=1:R
    for j=1:C
        if ass1(i,j)>0
            ass1(i,j)=1;
        end
        if ass2(i,j)>0
            ass2(i,j)=2;
        end
        if ass3(i,j)>0
            ass3(i,j)=3;
        end
    end
end

%representing the clustered image
fin1cluste=(ass1+ass2+ass3); %sum up the 3 labelled cluster
fin1cluste1=label2rgb(fin1cluste,@lines,'g','noshuffle');  %final segmented image   %%%%
axes(handles.axes4);
imshow(fin1cluste1)

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Image
I = im2double(Image);
I = imadjust(I,[.2 .3 0; .6 .7 1],[]);
I = imadjust(rgb2gray(I));
data = [I(:)]; % data array
[center,U,obj_fcn] = fcm(data,3,3); % Fuzzy C-means classification with 5 classes

% Finding the pixels for each class
maxU = max(U);
index1 = find(U(1,:) == maxU);
index2 = find(U(2,:) == maxU);
index3 = find(U(3,:) == maxU);

% Assigning pixel to each class by giving them a specific value
fcmImage(1:length(data))=0;       
fcmImage(index1)= 1;
fcmImage(index2)= 0.5;
fcmImage(index3)= 0;

% Reshapeing the array to a image
imagNew = reshape(fcmImage,size(I));
img = im2bw(imagNew,0);
axes(handles.axes5);
imshow(img);
