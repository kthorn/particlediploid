function varargout = particletrack2(varargin)
% PARTICLETRACK2 M-file for particletrack.fig
%      PARTICLETRACK2, by itself, creates a new PARTICLETRACK2 or raises the existing
%      singleton*.
%
%      H = PARTICLETRACK2 returns the handle to a new PARTICLETRACK2 or the handle to
%      the existing singleton*.
%
%      PARTICLETRACK2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARTICLETRACK2.M with the given input arguments.
%
%      PARTICLETRACK2('Property','Value',...) creates a new PARTICLETRACK2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before particletrack2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to particletrack2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help particletrack2

% Last Modified by GUIDE v2.5 05-May-2010 14:03:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @particletrack2_OpeningFcn, ...
    'gui_OutputFcn',  @particletrack2_OutputFcn, ...
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


% --- Executes just before particletrack2 is made visible.
function particletrack2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to particletrack2 (see VARARGIN)

% Choose default command line output for particletrack2
handles.output = hObject;
set(handles.remove_dots,'Enable','off');
set(handles.queue,'Enable','off');
set(handles.analyze,'Enable','off');

%default parameters
%both increased by two from haploid code
handles.data.len = 11; %size of box in XY to fit around each dot
handles.min_dot_dist = 13; %dots closer than this belong to same cell

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes particletrack2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = particletrack2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function dir_to_read_Callback(hObject, eventdata, handles)
% hObject    handle to dir_to_read (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get structure of data set
handles.data.MMstructure = MMstruct(get(handles.dir_to_read,'String'));
set(handles.wavelength_list,'String',handles.data.MMstructure.wavelengthlist);
guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of dir_to_read as text
%        str2double(get(hObject,'String')) returns contents of dir_to_read as a double


% --- Executes during object creation, after setting all properties.
function dir_to_read_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dir_to_read (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold as text
%        str2double(get(hObject,'String')) returns contents of threshold as a double


% --- Executes during object creation, after setting all properties.
function threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in analyze.
function analyze_Callback(hObject, eventdata, handles)
% hObject    handle to analyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nqueued = size(handles.data.queued,2);
for n = 1:nqueued
    set(handles.status,'String',['Analyzing ', sprintf('%d',n), ' of ', sprintf('%d',nqueued)]);
    drawnow expose update
    ims = MMparse(handles.data.queued(n).dir,[],{handles.data.queued(n).wavelength});
    model = fit_cell_models(ims, handles.data.queued(n).startcoords, handles.data.len);
    save(fullfile(handles.data.queued(n).dir,'analysis.mat'),'model','ims');
end
set(handles.status,'String','Done Analyzing');

guidata(hObject, handles);


% --- Executes on button press in queue.
function queue_Callback(hObject, eventdata, handles)
% hObject    handle to queue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles.data, 'queued')
    nqueued = size(handles.data.queued,2);
else
    nqueued = 0;
end
handles.data.queued(nqueued+1).startcoords = handles.data.startcoords;
handles.data.queued(nqueued+1).dir = get(handles.dir_to_read,'String');

wavelengthlist = get(handles.wavelength_list,'String');
wavelength = wavelengthlist{get(handles.wavelength_list,'Value')};
handles.data.queued(nqueued+1).wavelength = wavelength;

set(handles.status,'String','Dataset queued');

set(handles.analyze,'Enable','on');
guidata(hObject, handles);

% --- Executes on button press in find_dots.
function find_dots_Callback(hObject, eventdata, handles)
% hObject    handle to find_dots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%parameters
h=fspecial('log',5,1);
%remove points closer to edge of image stack than this
XYbound = 20;
Zbound = 4;

%image
t1=handles.data.images;
thresh = str2double(get(handles.threshold,'String'));

%filter image
%sharpen image to find dots
testfilt = sharpen_image(t1);

%find peaks
s=regionprops(testfilt>thresh,t1,'WeightedCentroid','MeanIntensity');
startcoords=zeros([size(s,1) 3]);
for n=1:size(s,1)
    startcoords(n,:)=round(s(n).WeightedCentroid);
end

%remove dots too close to edges of image or to other dots
d=squareform(pdist(startcoords));
[r,c]=find(d<handles.min_dot_dist & d>0);
pair_idx = find (r<c);
c=c(pair_idx);  %upper halves of close dots
r=r(pair_idx);

%iterate over paired dots, keeping brightest
del_list=[];
for n=1:size(r)
    if s(r(n)).MeanIntensity > s(c(n)).MeanIntensity
        del_list=[del_list; c(n)]; %r is brighter, so delete c
    else
        del_list=[del_list; r(n)];
    end
end

startcoords(del_list,:)=[]; %remove paired dots

XYbound=handles.data.len*2;
badcoords = any(startcoords(:,1:2)' <XYbound) | startcoords(:,1)'>512-XYbound | startcoords(:,2)'>512-XYbound;
startcoords(badcoords',:)=[];

handles.data.startcoords = startcoords;

update_image(handles.axes, max(t1,[],3), startcoords);
set(handles.status,'String',['Found ', sprintf('%d', size(startcoords,1)), ' Dots']);
set(handles.remove_dots,'Enable','on');
set(handles.queue,'Enable','on');
handles.data.time=1;

guidata(hObject, handles);


% --- Executes on button press in add_dots.
function add_dots_Callback(hObject, eventdata, handles)
% hObject    handle to add_dots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[x,y]=ginput(1);
x=round(x);
y=round(y);
[~, z]=max(handles.data.images(y,x,:,1));
handles.data.startcoords = [handles.data.startcoords; x y z];

update_image(handles.axes,max(handles.data.images(:,:,:,1),[],3),handles.data.startcoords);

set(handles.remove_dots,'Enable','on');
set(handles.queue,'Enable','on');

guidata(hObject, handles);


% --- Executes on button press in remove_dots.
function remove_dots_Callback(hObject, eventdata, handles)
% hObject    handle to remove_dots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[x,y]=ginput(1);
dist=zeros([1 size(handles.data.startcoords,1)]);
for n=1:size(handles.data.startcoords,1)
    dist(n) = (handles.data.startcoords(n,1)-x).^2 + (handles.data.startcoords(n,2)-y).^2;
end
[~, point_to_remove]=min(dist);
handles.data.startcoords(point_to_remove,:)=[];
update_image(handles.axes,max(handles.data.images(:,:,:,1),[],3),handles.data.startcoords);

guidata(hObject, handles);


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.status,'String','Loading Images');
drawnow expose update
wavelengthlist = get(handles.wavelength_list,'String');
wavelength = wavelengthlist{get(handles.wavelength_list,'Value')};

%load stack
handles.data.images=MMparse(get(handles.dir_to_read,'String'),1,{wavelength}); %load first time point

%background subtract
handles.data.images=handles.data.images-median(handles.data.images(:));

set(handles.status,'String',['Loaded ', sprintf('%d', size(handles.data.images,3) * size(handles.data.images,4)), ' Images']);
%clear values from previous run, if any
if isfield(handles.data,'model')
    handles.data = rmfield(handles.data,'model');
end
set(handles.queue,'Enable','off');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on selection change in wavelength_list.
function wavelength_list_Callback(hObject, eventdata, handles)
% hObject    handle to wavelength_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns wavelength_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from wavelength_list


% --- Executes during object creation, after setting all properties.
function wavelength_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to wavelength_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_image(target_axis,image,coords)
cla(target_axis)
imshow(image,[],'Parent',target_axis);
hold on
for n=1:size(coords,1)
    plot(coords(n,1),coords(n,2),'ro');
end
