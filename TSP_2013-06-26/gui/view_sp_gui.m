function varargout = view_sp_gui(varargin)
% VIEW_SP_GUI MATLAB code for view_sp_gui.fig
%      VIEW_SP_GUI, by itself, creates a new VIEW_SP_GUI or raises the existing
%      singleton*.
%
%      H = VIEW_SP_GUI returns the handle to a new VIEW_SP_GUI or the handle to
%      the existing singleton*.
%
%      VIEW_SP_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW_SP_GUI.M with the given input arguments.
%
%      VIEW_SP_GUI('Property','Value',...) creates a new VIEW_SP_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before view_sp_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to view_sp_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help view_sp_gui

% Last Modified by GUIDE v2.5 15-Apr-2013 11:34:25

global SPs f;
SPs = [];
f = 1;

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @view_sp_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @view_sp_gui_OutputFcn, ...
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

% --- Executes just before view_sp_gui is made visible.
function view_sp_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to view_sp_gui (see VARARGIN)

% Choose default command line output for view_sp_gui
handles.output = hObject;

handles = guidata(hObject);
handles.sp_labels = varargin{1};
handles.root = varargin{2};

if (numel(varargin)<3 || isempty(varargin{3}))
    folder = handles.root;
    handles.files = dir_images(folder);
else
    handles.files = varargin{3};
end

if (numel(varargin)<4 || isempty(varargin{4}))
    handles.frames = 1:numel(handles.files);
else
    handles.frames = varargin{4};
end

if (numel(varargin)<5 || isempty(varargin{5}))
	handles.flow = cell(size(handles.sp_labels,3));
    for i=1:size(handles.sp_labels,3)
        handles.flow{i} = zeros(size(handles.sp_labels,1),size(handles.sp_labels,2));
    end
else
    handles.flow = varargin{5};
end

if (numel(varargin)<6 || isempty(varargin{6}))
    [~,~,temp] = unique(handles.sp_labels);
else
    temp = varargin{6};
end

handles.sp_labels = reshape(temp, size(handles.sp_labels));
if (numel(varargin)<7 || isempty(varargin{7}))
    handles.sp_colors = distinguishable_colors(max(handles.sp_labels(:)));
else
    handles.sp_colors = varargin{7};
end



handles.f = 1;
handles.oim = im2double(imread([handles.root handles.files(handles.frames(handles.f)).name]));

% handles.flow = cell(1, size(sp_u,3));
% for f=1:size(sp_u,3)-1
%     handles.flow{f} = flowToColor(cat(3,handles.sp_u(:,:,f+1), handles.sp_u(:,:,f+1)));
% end
% handles.flow{end} = zeros(size(handles.oim));

handles.C = 50;
handles.SPs = [];
handles.SPsc = cell(handles.C,1);
handles.c = 1;
handles.colors = distinguishable_colors(50);
set(handles.txtColor, 'BackgroundColor', handles.colors(handles.c,:));

set(handles.slider_img, 'Value', 1);
set(handles.slider_img, 'Max', numel(handles.frames));
set(handles.slider_img, 'Min', 1);
set(handles.slider_img, 'SliderStep', [1/numel(handles.frames) 1/numel(handles.frames)]);


% Update handles structure
guidata(hObject, handles);

%guidata(hObject, guiData);
axes(handles.axes1);
hold off;
imagesc(handles.oim);
axis off;


if (false && ~isempty(handles.flow))
axes(handles.axes2);
hold off;
imagesc(handles.flow{handles.f});
axis off;
else
    temp_im = handles.sp_colors(handles.sp_labels(:,:,handles.f),:);
    temp_im = reshape(temp_im, size(handles.sp_labels,1),size(handles.sp_labels,2),3);
    imagesc(temp_im, 'Parent', handles.axes2);
end


% UIWAIT makes view_sp_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = view_sp_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Get default command line output from handles structure
%varargout{1} = handles.output;

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)








% --- Executes on button press in btn_single.
function btn_single_Callback(hObject, eventdata, handles)
% hObject    handle to btn_single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
display_borders_all(handles);
button = 1;

xdim = size(handles.sp_labels,1);
ydim = size(handles.sp_labels,2);

while (button==1)
    [y, x, button] = ginput(1);
    
    x = round(x);
    y = round(y);
    if (x>0 && x<=xdim && y>0 && y<=ydim)
        thislabel = handles.sp_labels(x, y,handles.f);
        if (button==1)
            if (any(handles.SPs==thislabel))
                for i=1:handles.C
                    handles.SPsc{i}(handles.SPsc{i}==thislabel) = [];
                end
            else
                handles.SPs(end+1) = thislabel;
            end
            handles.SPsc{handles.c}(end+1) = thislabel;
            display_borders_all(handles);
        end
    end
end
display_borders_cur(handles);
guidata(hObject, handles);



















% --- Executes on button press in btn_line.
function btn_line_Callback(hObject, eventdata, handles)
% hObject    handle to btn_line (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
display_borders_all(handles);
[points] = imlinesegment();

xdim = size(handles.sp_labels,1);
ydim = size(handles.sp_labels,2);

N = size(points,1);
for i=1:N-1
    xdelta = abs(points(i,2) - points(i+1,2));
    ydelta = abs(points(i,1) - points(i+1,1));

    nmid = round(max(ydelta,xdelta));
    x = round(linspace(points(i,2),points(i+1,2),nmid));
    y = round(linspace(points(i,1),points(i+1,1),nmid));
    
    mask = x>0 & x<=xdim & y>0 & y<=ydim;
    x = x(mask);
    y = y(mask);

    values = unique(handles.sp_labels(sub2ind([xdim,ydim,handles.f], x, y, repmat(handles.f,size(x)))));
    for thislabel=values
        if (any(handles.SPs==thislabel))
            for c=1:handles.C
                handles.SPsc{c}(handles.SPsc{c}==thislabel) = [];
            end
        else
            handles.SPs(end+1) = thislabel;
        end
        handles.SPsc{handles.c}(end+1) = thislabel;
    end
    %handles.SPs(end+1:end+numel(values)) = values;
end
%handles.SPs = unique(handles.SPs);
display_borders_cur(handles);
guidata(hObject, handles);


            
           



% --- Executes on button press in btn_polygon.
function btn_polygon_Callback(hObject, eventdata, handles)
% hObject    handle to btn_polygon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
display_borders_all(handles);
h = impoly();
mask = createMask(h);
values = unique(handles.sp_labels(find(mask)+(handles.f-1)*size(handles.oim,1)*size(handles.oim,2)))';
for thislabel=values
    if (any(handles.SPs==thislabel))
        for c=1:handles.C
            handles.SPsc{c}(handles.SPsc{c}==thislabel) = [];
        end
    else
        handles.SPs(end+1) = thislabel;
    end
    handles.SPsc{handles.c}(end+1) = thislabel;
end
%handles.SPs(end+1:end+numel(values)) = values;
delete(h);
%handles.SPs = unique(handles.SPs);
display_borders_cur(handles);
guidata(hObject, handles);


% --- Executes on button press in btn_rectangle.
function btn_rectangle_Callback(hObject, eventdata, handles)
% hObject    handle to btn_rectangle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
display_borders_all(handles);
h = imrect();
mask = createMask(h);
values = unique(handles.sp_labels(find(mask)+(handles.f-1)*size(handles.oim,1)*size(handles.oim,2)))';
for thislabel=values
    if (any(handles.SPs==thislabel))
        for c=1:handles.C
            handles.SPsc{c}(handles.SPsc{c}==thislabel) = [];
        end
    else
        handles.SPs(end+1) = thislabel;
    end
    handles.SPsc{handles.c}(end+1) = thislabel;
end
%handles.SPs(end+1:end+numel(values)) = values;
delete(h);
%handles.SPs = unique(handles.SPs);
display_borders_cur(handles);
guidata(hObject, handles);


% --- Executes on button press in btn_freehand.
function btn_freehand_Callback(hObject, eventdata, handles)
% hObject    handle to btn_freehand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
display_borders_all(handles);
h = imfreehand();
mask = createMask(h);
values = unique(handles.sp_labels(find(mask)+(handles.f-1)*size(handles.oim,1)*size(handles.oim,2)))';
for thislabel=values
    if (any(handles.SPs==thislabel))
        for c=1:handles.C
            handles.SPsc{c}(handles.SPsc{c}==thislabel) = [];
        end
    else
        handles.SPs(end+1) = thislabel;
    end
    handles.SPsc{handles.c}(end+1) = thislabel;
end
%handles.SPs(end+1:end+numel(values)) = values;
delete(h);
%handles.SPs = unique(handles.SPs);
display_borders_cur(handles);
guidata(hObject, handles);


% --- Executes on button press in btn_delete.
function btn_delete_Callback(hObject, eventdata, handles)
% hObject    handle to btn_delete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
display_borders_all(handles);
button = 1;
while (button==1)
    [y, x, button] = ginput(1);
    if (button==1 && any(handles.SPs(:)==handles.sp_labels(round(x), round(y),handles.f)))
        thislabel = handles.sp_labels(round(x), round(y), handles.f);
        if (any(handles.SPs==thislabel))
            handles.SPs(handles.SPs==thislabel) = [];
            for c=1:handles.C
                handles.SPsc{c}(handles.SPsc{c}==thislabel) = [];
            end
        end
        display_borders_all(handles);
    end
end
display_borders_cur(handles);
guidata(hObject, handles);


% --- Executes on button press in btn_deleteall.
function btn_deleteall_Callback(hObject, eventdata, handles)
% hObject    handle to btn_deleteall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SPs = [];
handles.SPsc = cell(handles.C,1);
guidata(hObject, handles);
display_borders_cur(handles);



function btnTransparent_Callback(hObject,eventdata, handles)
display_borders_cur(handles);

function btnFilled_Callback(hObject,eventdata, handles)
display_borders_cur(handles);

function btnColorCoded_Callback(hObject,eventdata, handles)
display_borders_cur(handles);


function [im] = display_borders_all(handles)
mask = ismember(handles.sp_labels(:,:,handles.f), handles.SPs);
im = setPixelColors(handles.oim, find(mask), [0 1 0]);

borders = is_border_valsIMPORT(double(handles.sp_labels(:,:,handles.f)));
im = setPixelColors(im, find(borders), [1 0 0]);
axes(handles.axes1);
hold off;
image(im, 'Parent', handles.axes1);
axis off;

function [im] = display_borders_cur(handles)
im = handles.oim;

im2 = handles.sp_colors(handles.sp_labels(:,:,handles.f),:);
im2 = reshape(im2, size(handles.sp_labels,1),size(handles.sp_labels,2),3);

N = size(im,1)*size(im,2);
transparent = 1-get(handles.slider_transparency,'Value');

bgmask = ~ismember(handles.sp_labels(:,:,handles.f), handles.SPs);
for c=1:handles.C
    if (~isempty(handles.SPsc{c}))
        color = handles.colors(c,:);
        cmask = ismember(handles.sp_labels(:,:,handles.f), handles.SPsc{c});
        indices = find(cmask);
        im(indices) = im(indices)*transparent + (1-transparent)*color(1);
        im(indices+N) = im(indices+N)*transparent + (1-transparent)*color(2);
        im(indices+2*N) = im(indices+2*N)*transparent + (1-transparent)*color(3);

        
        borders = is_border_valsIMPORT(double(handles.sp_labels(:,:,handles.f)).*cmask);
        borders = borders & (cmask | bgmask);
        im = setPixelColors(im, find(borders), color);
        im2 = setPixelColors(im2, find(cmask), color);
        im2 = setPixelColors(im2, find(borders), [1 1 1]);
    end
end
%axes(handles.axes1);
%hold off;
image(im, 'Parent', handles.axes1);
%axis off;

if (false && ~isempty(handles.flow))
%axes(handles.axes2);
%hold off;
imagesc(handles.flow{handles.f}, 'Parent', handles.axes2);
%axis off;
else
    temp_im = handles.sp_colors(handles.sp_labels(:,:,handles.f),:);
    temp_im = reshape(temp_im, size(handles.sp_labels,1),size(handles.sp_labels,2),3);
    imagesc(im2, 'Parent', handles.axes2);
end



% --- Executes on button press in btn_next.
function btn_next_Callback(hObject, eventdata, handles)
% hObject    handle to btn_next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.f < numel(handles.frames))
    handles.f = handles.f + 1;
    handles.oim = im2double(imread([handles.root handles.files(handles.frames(handles.f)).name]));
    display_borders_cur(handles);
    set(handles.slider_img, 'Value', handles.f);
    guidata(hObject, handles);
end


% --- Executes on button press in btn_previous.
function btn_previous_Callback(hObject, eventdata, handles)
% hObject    handle to btn_previous (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (handles.f>1)
    handles.f = handles.f - 1;
    handles.oim = im2double(imread([handles.root handles.files(handles.frames(handles.f)).name]));
    display_borders_cur(handles);
    set(handles.slider_img, 'Value', handles.f);
    guidata(hObject, handles);
end


% --- Executes on slider movement.
function slider_img_Callback(hObject, eventdata, handles)
% hObject    handle to slider_img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.f = round(get(hObject,'Value'));
set(hObject, 'Value', handles.f);

handles.oim = im2double(imread([handles.root handles.files(handles.frames(handles.f)).name]));
display_borders_cur(handles);
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function slider_img_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_img (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on scroll wheel click while the figure is in focus.
function figure1_WindowScrollWheelFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	VerticalScrollCount: signed integer indicating direction and number of clicks
%	VerticalScrollAmount: number of lines scrolled for each click
% handles    structure with handles and user data (see GUIDATA)
if (handles.f < numel(handles.frames) && eventdata.VerticalScrollCount==1)
    handles.f = handles.f + 1;
    handles.oim = im2double(imread([handles.root handles.files(handles.frames(handles.f)).name]));
    display_borders_cur(handles);
    set(handles.slider_img, 'Value', handles.f);
    guidata(hObject, handles);
elseif (handles.f>1 && eventdata.VerticalScrollCount==-1)
    handles.f = handles.f - 1;
    handles.oim = im2double(imread([handles.root handles.files(handles.frames(handles.f)).name]));
    display_borders_cur(handles);
    set(handles.slider_img, 'Value', handles.f);
    guidata(hObject, handles);
end


% --- Executes on slider movement.
function slider_transparency_Callback(hObject, eventdata, handles)
% hObject    handle to slider_transparency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
display_borders_cur(handles);


% --- Executes during object creation, after setting all properties.
function slider_transparency_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_transparency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function SaveMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

count = 1;
for f=handles.frames
    disp(['Saving file ' num2str(count) ' / ' num2str(numel(handles.frames))]);
    
    
    im = im2double(imread([handles.root handles.files(f).name]));
    N = size(im,1)*size(im,2);
    transparent = 1-get(handles.slider_transparency,'Value');

    bgmask = ~ismember(handles.sp_labels(:,:,f), handles.SPs);
    for c=1:handles.C
        if (~isempty(handles.SPsc{c}))
            color = handles.colors(c,:);
            cmask = ismember(handles.sp_labels(:,:,f), handles.SPsc{c});
            indices = find(cmask);
            im(indices) = im(indices)*transparent + (1-transparent)*color(1);
            im(indices+N) = im(indices+N)*transparent + (1-transparent)*color(2);
            im(indices+2*N) = im(indices+2*N)*transparent + (1-transparent)*color(3);


            borders = is_border_valsIMPORT(double(handles.sp_labels(:,:,f)).*cmask);
            borders = borders & (cmask | bgmask);
            im = setPixelColors(im, find(borders), color);
        end
    end
    
    imwrite(im2uint8(im), ['saves/' handles.files(f).name(1:end-4) '.png']);
    count = count + 1;
end


% --- Executes on button press in btnMinus.
function btnMinus_Callback(hObject, eventdata, handles)
% hObject    handle to btnMinus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.c = max(1,handles.c - 1);
set(handles.txtColor, 'BackgroundColor', handles.colors(handles.c,:));
guidata(hObject, handles);


% --- Executes on button press in btnPlus.
function btnPlus_Callback(hObject, eventdata, handles)
% hObject    handle to btnPlus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.c = min(handles.C,handles.c + 1);
set(handles.txtColor, 'BackgroundColor', handles.colors(handles.c,:));
guidata(hObject, handles);


% --- Executes on button press in btn_all.
function btn_all_Callback(hObject, eventdata, handles)
% hObject    handle to btn_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

values = unique(handles.sp_labels)';
for c=1:handles.C
    handles.SPsc{c} = [];
end
handles.SPsc{handles.c} = values;
handles.SPs = values;
display_borders_cur(handles);
guidata(hObject, handles);


% --- Executes on button press in btn_frame.
function btn_frame_Callback(hObject, eventdata, handles)
% hObject    handle to btn_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

values = unique(handles.sp_labels(:,:,handles.f))';
for thislabel=values
    if (any(handles.SPs==thislabel))
        for c=1:handles.C
            handles.SPsc{c}(handles.SPsc{c}==thislabel) = [];
        end
    else
        handles.SPs(end+1) = thislabel;
    end
    handles.SPsc{handles.c}(end+1) = thislabel;
end
display_borders_cur(handles);
guidata(hObject, handles);
