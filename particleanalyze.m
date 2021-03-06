function varargout = particleanalyze(varargin)
% PARTICLEANALYZE M-file for particleanalyze.fig
%      PARTICLEANALYZE, by itself, creates a new PARTICLEANALYZE or raises the existing
%      singleton*.
%
%      H = PARTICLEANALYZE returns the handle to a new PARTICLEANALYZE or the handle to
%      the existing singleton*.
%
%      PARTICLEANALYZE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PARTICLEANALYZE.M with the given input arguments.
%
%      PARTICLEANALYZE('Property','Value',...) creates a new PARTICLEANALYZE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before particleanalyze_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to particleanalyze_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help particleanalyze

% Last Modified by GUIDE v2.5 05-Nov-2010 15:57:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @particleanalyze_OpeningFcn, ...
    'gui_OutputFcn',  @particleanalyze_OutputFcn, ...
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


% --- Executes just before particleanalyze is made visible.
function particleanalyze_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to particleanalyze (see VARARGIN)

% Choose default command line output for particleanalyze
handles.output = hObject;

disp('Initializing ... ')

handles.data.model=evalin('base','model');

if evalin('base', 'exist(''ims'', ''var'')');
    handles.data.ims=evalin('base','ims');
else    
    dir = uigetdir('Z:\','Image directory?');
    handles.data.ims = MMparse(dir);
end
handles.data.time=1;
handles.data.selected=-1;

disp('Modeling intensities ... ')
%calculate intensities
%model each dot independently

for n = 1:size(handles.data.model,2)
    handles.data.model(n).flags={};
    int_2dot = handles.data.model(n).params(:,13)+handles.data.model(n).params(:,17);
    [d1,d2] = reassign_dots(handles.data.model(n).params);
    handles.data.model(n).dot1_I = d1(:,4);
    handles.data.model(n).dot2_I = d2(:,4);
    
    %model intensity disappearance
    handles.data.model(n).model1_I = fit_disappearance(handles.data.model(n).dot1_I');
    handles.data.model(n).model2_I = fit_disappearance(handles.data.model(n).dot2_I');
    
    %find potential disappearing dots
    if (handles.data.model(n).model1_I(4)/handles.data.model(n).model1_I(1) < 0.3)...
            && (handles.data.model(n).model2_I(4)/handles.data.model(n).model2_I(1) < 0.3)
        handles.data.model = addflag(handles.data.model, n, 'curious');
    end
end

update_image(handles.axes,handles.data);

n_cells = size(handles.data.model,2);
set(handles.cell_slider,'Max',n_cells);
set(handles.cell_slider,'Min',1);
set(handles.cell_slider,'Value',1);

set(handles.cell_slider,'SliderStep',[1/n_cells 10/n_cells]);
disp('done')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes particleanalyze wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = particleanalyze_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in plot_intensity.
function plot_intensity_Callback(hObject, eventdata, handles)
% hObject    handle to plot_intensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.data.selected > 0
    int_1dot = handles.data.model(handles.data.selected).dot1_I;
    int_2dot = handles.data.model(handles.data.selected).dot2_I;    
    figure(1)
    clf
    subplot(2,1,1)
    plot(int_1dot);
    hold on
    plot(model_results(handles.data.model(handles.data.selected).model1_I),'r')
    title('1st dot intensity')
    subplot(2,1,2)
    plot(int_2dot);
    hold on
    plot(model_results(handles.data.model(handles.data.selected).model2_I),'r')
    title('2nd dot intensity')
end
guidata(hObject, handles);

% --- Executes on mouse press over axes background.
function axes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in export.
function export_Callback(hObject, eventdata, handles)
% hObject    handle to export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

results = handles.data.model;
results = rmfield(results, {'im', 'fitim', 'raw_params'});

variable_name=inputdlg('Name for output data');
assignin('base', variable_name{1}, results);


% --- Executes on button press in selectcell.
function selectcell_Callback(hObject, eventdata, handles)
% hObject    handle to selectcell (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[x,y]=ginput(1);
dist=zeros([1 size(handles.data.model,2)]);
for n=1:size(handles.data.model,2)
    dist(n) = (handles.data.model(n).params(1,10)-x).^2 + (handles.data.model(n).params(1,11)-y).^2;
end
[junk, selected_point]=min(dist);
handles = update_selection(handles, selected_point);

guidata(hObject, handles);


% --- Executes on slider movement.
function cell_slider_Callback(hObject, eventdata, handles)
% hObject    handle to cell_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

n_cell = round(get(hObject,'Value'));
set(handles.cell_slider, 'Value', n_cell);
set(handles.cell_text,'String',sprintf('%u',n_cell));
handles = update_selection(handles, n_cell);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cell_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function cell_text_Callback(hObject, eventdata, handles)
% hObject    handle to cell_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cell_text as text
%        str2double(get(hObject,'String')) returns contents of cell_text as a double

n_cell = str2double(get(handles.cell_text,'String'));
set(handles.cell_slider, 'Value', n_cell);
handles = update_selection(handles, n_cell);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function cell_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in showmodel.
function showmodel_Callback(hObject, eventdata, handles)
% hObject    handle to showmodel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.data.selected > 0
    %show one and two dot fits and relative performance    
    figure(1)
    clf
 
    for t=1:60
        subplot(4,15,t)
        dispim = pair_2dot_images(handles.data.model, handles.data.selected, t);
        imshow(dispim, [], 'InitialMagnification', 'fit')
        title(strcat(sprintf('%d',t)));
    end
end

    
% --- Executes on button press in play_tracks.
function play_tracks_Callback(hObject, eventdata, handles)
% hObject    handle to play_tracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for t=1:size(handles.data.ims,4)
    handles.data.time=t;
    update_image(handles.axes,handles.data);
    
    drawnow expose update
end


% --- Executes on button press in delete_tracks.
function delete_tracks_Callback(hObject, eventdata, handles)
% hObject    handle to delete_tracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.model(handles.data.selected)=[];
update_image(handles.axes,handles.data);

guidata(hObject, handles);


% --- Executes on button press in curious_flag.
function curious_flag_Callback(hObject, eventdata, handles)
% hObject    handle to curious_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of curious_flag

if get(hObject, 'Value') == get(hObject, 'Max')
    handles.data.model = addflag(handles.data.model, handles.data.selected, 'curious');
elseif get(hObject, 'Value') == get(hObject, 'Min')    
    handles.data.model = remflag(handles.data.model, handles.data.selected, 'curious');
end

guidata(hObject, handles);



% --- Executes on button press in disappearing_flag.
function disappearing_flag_Callback(hObject, eventdata, handles)
% hObject    handle to disappearing_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of disappearing_flag

if get(hObject, 'Value') == get(hObject, 'Max')
    handles.data.model = addflag(handles.data.model, handles.data.selected, 'disappearing');
elseif get(hObject, 'Value') == get(hObject, 'Min')    
    handles.data.model = remflag(handles.data.model, handles.data.selected, 'disappearing');
end

guidata(hObject, handles);


% --- Executes on button press in single_flag.
function single_flag_Callback(hObject, eventdata, handles)
% hObject    handle to single_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of single_flag
if get(hObject, 'Value') == get(hObject, 'Max')
    handles.data.model = addflag(handles.data.model, handles.data.selected, 'single');
elseif get(hObject, 'Value') == get(hObject, 'Min')    
    handles.data.model = remflag(handles.data.model, handles.data.selected, 'single');
end

guidata(hObject, handles);


% --- Executes on button press in splitting_flag.
function splitting_flag_Callback(hObject, eventdata, handles)
% hObject    handle to splitting_flag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of splitting_flag
if get(hObject, 'Value') == get(hObject, 'Max')
    handles.data.model = addflag(handles.data.model, handles.data.selected, 'splitting');
elseif get(hObject, 'Value') == get(hObject, 'Min')    
    handles.data.model = remflag(handles.data.model, handles.data.selected, 'splitting');
end

guidata(hObject, handles);

function crossover_Callback(hObject, eventdata, handles)
% hObject    handle to crossover (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of crossover as text
%        str2double(get(hObject,'String')) returns contents of crossover as a double


% --- Executes during object creation, after setting all properties.
function crossover_CreateFcn(hObject, eventdata, handles)
% hObject    handle to crossover (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function penalty_Callback(hObject, eventdata, handles)
% hObject    handle to penalty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of penalty as text
%        str2double(get(hObject,'String')) returns contents of penalty as a double


% --- Executes during object creation, after setting all properties.
function penalty_CreateFcn(hObject, eventdata, handles)
% hObject    handle to penalty (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in recall_dots.
function recall_dots_Callback(hObject, eventdata, handles)
% hObject    handle to recall_dots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Determining number of dots ... ')
crossover = str2double(get(handles.crossover, 'String'));
penalty = str2double(get(handles.penalty, 'String'));
for n=1:size(handles.data.model,2)
    handles.data.model(n).ndots = 2;   
end
update_image(handles.axes,handles.data);
disp('Done')
guidata(hObject, handles);


function update_image(target_axis,data_struct)
time=data_struct.time;
cla(target_axis)
imshow(max(data_struct.ims(:,:,:,time),[],3),[],'Parent',target_axis);
hold on
for n=1:size(data_struct.model,2)
    color = 'r';
    if any (strcmp(data_struct.model(n).flags, 'curious'))
        color = 'g';
    end
    if (n == data_struct.selected)
        color = 'w';
    end
    plot(target_axis,data_struct.model(n).params(time,10),data_struct.model(n).params(time,11),[color,'s']);
    plot(target_axis,data_struct.model(n).params(time,14),data_struct.model(n).params(time,15),[color,'o']);
end

function dispim = pair_2dot_images(model, n, t)
dispim=max(model(n).im(:,:,:,t),[],3);
dispim=dispim-min(dispim(:));
dispim=[dispim;zeros([1 size(dispim,1)])];
fitim=max(model(n).fitim(:,:,:,t),[],3);
fitim=fitim-min(fitim(:));
dispim=[dispim; fitim];

function model = addflag (model, n, flag)
model(n).flags = [model(n).flags, flag];

function model = remflag (model, n, flag)
model(n).flags(strcmp(model(n).flags, flag)) = [];

function handles = update_selection(handles, selected_cell)

handles.data.selected=selected_cell;
set(handles.cell_text, 'String', sprintf('%u', selected_cell));
set(handles.cell_slider, 'Value', selected_cell);
if any(strcmp(handles.data.model(selected_cell).flags, 'curious'))
    set(handles.curious_flag, 'Value', get(handles.curious_flag, 'Max'));
else    
    set(handles.curious_flag, 'Value', get(handles.curious_flag, 'Min'));
end
if any(strcmp(handles.data.model(selected_cell).flags, 'single'))
    set(handles.single_flag, 'Value', get(handles.single_flag, 'Max'));
else    
    set(handles.single_flag, 'Value', get(handles.single_flag, 'Min'));
end
if any(strcmp(handles.data.model(selected_cell).flags, 'splitting'))
    set(handles.splitting_flag, 'Value', get(handles.splitting_flag, 'Max'));
else
    set(handles.splitting_flag, 'Value', get(handles.splitting_flag, 'Min'));
end
if any(strcmp(handles.data.model(selected_cell).flags, 'disappearing'))
    set(handles.disappearing_flag, 'Value', get(handles.disappearing_flag, 'Max'));
else
    set(handles.disappearing_flag, 'Value', get(handles.disappearing_flag, 'Min'));
end
update_image(handles.axes,handles.data);

function [new_dot1, new_dot2] = reassign_dots(params)
%returns new_dot1, new_dot2 which are dot1 and dot2 from the original
%params sorted to minimize the distance change between time points

prev_dot1 = params(1,10:12);
prev_dot2 = params(1,14:16);
new_dot1(1,:) = params(1,10:13);
new_dot2(1,:) = params(1,14:17);

for t = 2:size(params,1)
    dot1 = params(t,10:12);
    dot2 = params(t,14:16);
    dist11 = sqrt(sum((prev_dot1 - dot1).^2));
    dist12 = sqrt(sum((prev_dot1 - dot2).^2));
    dist21 = sqrt(sum((prev_dot2 - dot1).^2));
    dist22 = sqrt(sum((prev_dot2 - dot2).^2));
    
    if (dist11+dist22) < (dist12+dist21);
        new_dot1(t,:) = params(t,10:13);
        new_dot2(t,:) = params(t,14:17);
    else
        new_dot1(t,:) = params(t,14:17);
        new_dot2(t,:) = params(t,10:13);
    end
    prev_dot1 = new_dot1(t,1:3);
    prev_dot2 = new_dot2(t,1:3);
end
