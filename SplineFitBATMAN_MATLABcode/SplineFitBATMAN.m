function varargout = SplineFitBATMAN(varargin)
% written 120213 Dr. Jie Hao, Imperial College London

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @SplineFitBATMAN_OpeningFcn, ...
    'gui_OutputFcn',  @SplineFitBATMAN_OutputFcn, ...
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


% --- Executes just before SplineFitBATMAN is made visible.
function SplineFitBATMAN_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for SplineFitBATMAN
handles.output = hObject;
handles.multi = cell(0);
handles.multi_ID = 0;
handles.nmrdata = [];
handles.ppm = [];
handles.dorder = [];
handles.flagpoint = 3;
handles.mksize = 12;
handles.ppmx = [];
handles.plotFigure = [];
handles.axHdl = [];
handles.x = [];
handles.y = [];
handles.h2 = [];
handles.hspline = [];
handles.hselect = [];
handles.a = [];
% axes(handles.logo);
% Update handles structure
guidata(hObject, handles);


% UIWAIT makes SplineFitBATMAN wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SplineFitBATMAN_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function BrowseNMR_Callback(hObject, eventdata, handles)
[FileName,PathName,FilterIndex]  = uigetfile('*.txt', 'Select NMR spectra file');
set(handles.NMRdir,'String',fullfile(PathName, FileName));


function NMRdir_Callback(hObject, eventdata, handles)

function NMRdir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function chemshiftDir_Callback(hObject, eventdata, handles)

function chemshiftDir_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function BrowseChemShift_Callback(hObject, eventdata, handles)
[FileName,PathName,FilterIndex]  = uigetfile('/Users/jhao2/Desktop/runBATMAN/BatmanInput/*.csv', 'Select chemical shift per spectrum file');
set(handles.chemshiftDir,'String',fullfile(PathName, FileName));

% function ppmZoomFrom_Callback(hObject, eventdata, handles)

function ppmZoomFrom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% function ppmZoomReset_Callback(hObject, eventdata, handles)
% replotzoomin(hObject, handles);

% function ppmZoomTo_Callback(hObject, eventdata, handles)

function ppmZoomTo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function multList_Callback(hObject, eventdata, handles)
tmp = get(hObject,'Value');
% if handles.multi_ID ~= tmp
handles.multi_ID = tmp;
% else
%     warndlg('Please choose a different multiplet, or click "Select and Run" to repeat with the same multiplet.','!! Warning !!');
% end
guidata(hObject, handles);

function multList_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minpoint_Callback(hObject, eventdata, handles)

function ppmRange_Callback(hObject, eventdata, handles)

function ppmRange_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function maxpoint_Callback(hObject, eventdata, handles)

function splinepoint_Callback(hObject, eventdata, handles)

function runSpline_Callback(hObject, eventdata, handles)
[hObject, handles] = deletePlot(hObject, handles);

if (handles.multi_ID == 1 || handles.multi_ID == 0)
    warndlg('Please choose a multiplet.','!! Warning !!');
    %     multList_Callback(hObject, eventdata, handles);
    return;
end
set(0,'CurrentFigure', handles.plotFigure);
% set(handles.plotFigure,'CurrentAxes', handles.axHdl);
[x,y,button] = ginput;
if (length(x) <2)
    warndlg('Please select more than 2 points. Click "Select and Run" to re-select.','!! Warning !!');
    [hObject, handles] = deletePlot(hObject, handles);
    return;
end
handles.x = x;
handles.y = y;

[hObject, handles] = replotspinepoints(hObject,handles);
guidata(hObject, handles);


function saveNext_Callback(hObject, eventdata, handles)
[hObject, handles] = deletePlot(hObject, handles);

m = handles.ppmx(1,:);
if length(m) == (size(handles.multi,2)-2)
    [~,d] = sort(handles.dorder);
    m = m(d);
    for i = 3:size(handles.multi,2)
        handles.multi{handles.multi_ID,i} = num2str(m(i-2));
    end
else
    warndlg('Number of spectra mismatch number of points selected. Please re-select. Click on points from (at least) the first and last spectra.');
    [hObject, handles] = deletePlot(hObject, handles);
    return;
end
guidata(hObject, handles);

function savewrite_Callback(hObject, eventdata, handles)
m = handles.ppmx(1,:);
if length(m) == (size(handles.multi,2)-2)
    [~,d] = sort(handles.dorder);
    m = m(d);
    for i = 3:size(handles.multi,2)
        handles.multi{handles.multi_ID,i} = num2str(m(i-2));
    end
else
    warndlg('Number of spectra mismatch number of points selected. Please re-select. Click on points from (at least) the first and last spectra.');
    [hObject, handles] = deletePlot(hObject, handles);
    return;
end
guidata(hObject, handles);
write_mixed_csv(get(handles.chemshiftDir,'String'), handles.multi);
msgbox('Chemical shift per spectrum saved in file.');


function ppmSortFrom_Callback(hObject, eventdata, handles)


function ppmSortFrom_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function yError_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function offset_Callback(hObject, eventdata, handles)

function offset_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Sortplotstack_Callback(hObject, eventdata, handles)
a = read_mixed_csv(get(handles.chemshiftDir,'String'),',');
for n = 2:size(a,1)
    b{1,n} =[a{n,1} ' ' a{n,2}];
end
set(handles.multList,'String',b);
handles.multi = a;

ppmSortFrom = str2num(get(handles.ppmSortFrom,'String'));
ppmSortTo = str2num(get(handles.ppmSortTo,'String'));

sortlimit = sort([ppmSortFrom, ppmSortTo]);

Data = readNMRdata(get(handles.NMRdir,'String'));
data = Data(:,2:end)';
if ((size(a,2)-2)~=size(data,1))
    warndlg('Number of spectra mismatch in .txt file and .csv file. Please check input files.');
    return;
end

ppm=Data(:,1)';

handles.nmrdata = data;
handles.ppm = ppm;
% handles.Xn = data(dorder,:);
ind = find(ppm<= sortlimit(2) & ppm>=sortlimit(1));
[~, bid]=max(data(:,ind)');
[~, dorder]= sort(bid);
dorder = fliplr(dorder);
handles.dorder = dorder;
offset = str2num(get(handles.offset,'String'));
[n m]=size(data);
% ph = (1:m)/10;
for i = 1:length(dorder)%size(handles.multi,2)
    lab{1,i} =handles.multi{1,dorder(i)+2};
end
if (~isempty(handles.plotFigure))
    set(0,'CurrentFigure', handles.plotFigure);
    handles.a = axis;
end
handles.plotFigure = figure(1);
clf
stackplot(data(dorder,:),ppm,lab,offset*ones(m,1),0);
set(gca,'XDir','reverse');
if (~(isempty(handles.a)))
    if (handles.a(1) ~= 0)
        xlim(handles.a(1:2));
        index = find(ppm<=handles.a(2) & ppm>=handles.a(1));
        ylimits=[-0.05 n*offset+max(data(dorder(end),index))];
        ylim(ylimits);
    end
end

guidata(hObject, handles);

function ppmSortTo_Callback(hObject, eventdata, handles)

function ppmSortTo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function uipanel9_SelectionChangeFcn(hObject, eventdata, handles)
if hObject == handles.maxpoint
    flagpoint = 1;
elseif hObject == handles.minpoint
    flagpoint = 2;
else
    flagpoint = 3;
end
handles.flagpoint = flagpoint;
guidata(hObject, handles);


function relocatep_Callback(hObject, eventdata, handles)
[hObject, handles] = deletePlot(hObject, handles);
[hObject, handles] = replotspinepoints(hObject, handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function logo_CreateFcn(hObject, eventdata, handles)
axes(hObject)
imshow('batmanlogo.png')
