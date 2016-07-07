function varargout = viewer(varargin)
% VIEWER MATLAB code for viewer.fig
%      VIEWER, by itself, creates a new VIEWER or raises the existing
%      singleton*.
%
%      H = VIEWER returns the handle to a new VIEWER or the handle to
%      the existing singleton*.
%
%      VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWER.M with the given input arguments.
%
%      VIEWER('Property','Value',...) creates a new VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before viewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help viewer

% Last Modified by GUIDE v2.5 16-Feb-2016 12:31:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @viewer_OpeningFcn, ...
                   'gui_OutputFcn',  @viewer_OutputFcn, ...
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


% --- Executes just before viewer is made visible.
function viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to viewer (see VARARGIN)
% if nargin ~= 6
%     disp('ERROR: expected object arrays EVENTMATRIX, MUTATIONS, STRAINS in this order');
% end

colors=load('mycolors');
handles.colors=colors;
handles.allparams = varargin{1};
% EVENTMATRIX=varargin{1};
% MUTATIONS=varargin{2};
% STRAINS=varargin{3};
% handles.EVENTMATRIX = EVENTMATRIX;
% handles.MUTATIONS = MUTATIONS;
% handles.STRAINS = STRAINS;

% initally just view a calls heatmap
blue=[55 67 153]/255;
lyel=[255 255 135]/255;
dyel=[255 191 0]/255;
cmap=flipud([lyel; blue; dyel; 0 0 0; .3 .3 .3]);
% xlab = {STRAINS.name};
% ylab = [MUTATIONS.gene_product];
% heatmap(param,xlab,ylab,[],'colormap',cmap,'showallticks',1,'tickangle',90);
% heatmap(handles.allparams.mutmat,{},{},[],'colormap',cmap,'showallticks',1,'tickangle',90);

listboxparams = {...
    'consensus',...
    'polymorphism',...
    'strandBias',...
    'ksTest',...
    'biasE',...
    'biasP',...
    'majorAFreq',...
    'covTotal',...
    'covTop',...
    'covBot',...
    'freq',...
    'consensusFC',...
    'polymorphismFC',...
    'strandBiasFC',...
    'ksTestFC',...
    'biasEFC',...
    'biasPFC',...
    'majorAFreqFC',...
    'covTotalFC',...
    'covTopFC',...
    'covBotFC',...
    'freqFC'};
handles.listboxparams = listboxparams;
    
handles.ParamsListX.String = listboxparams;
handles.ParamsListY.String = listboxparams;
% Choose default command line output for viewer
handles.output = hObject;
handles.axislabelX = '';
handles.axislabelY = '';
handles.xvals = [];
handles.yvals = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = viewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --- Executes during object creation, after setting all properties.
function ParamsListX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ParamsListX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function ParamsListY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ParamsListY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ParamsListX.
function ParamsListX_Callback(hObject, eventdata, handles)
% hObject    handle to ParamsListX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns ParamsListX contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ParamsListX
val = handles.ParamsListX.Value;
% set(handles, 'currentaxes', handles.MainPlot);
handles.figure1.CurrentAxes = handles.MainPlot;
handles.axislabelX = handles.listboxparams{val};
xvals = getvalue(hObject, eventdata, handles, val);
handles.xvals = xvals;
guidata(hObject, handles);
mainAxUpdate(hObject, eventdata, handles)

% --- Executes on selection change in ParamsListY.
function ParamsListY_Callback(hObject, eventdata, handles)
% hObject    handle to ParamsListY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val = handles.ParamsListY.Value;
% set(handles, 'currentaxes', handles.MainPlot);
handles.figure1.CurrentAxes = handles.MainPlot;
handles.axislabelY = handles.listboxparams{val};
% Hints: contents = cellstr(get(hObject,'String')) returns ParamsListY contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ParamsListY
yvals = getvalue(hObject, eventdata, handles, val);
handles.yvals = yvals;
guidata(hObject, handles);
mainAxUpdate(hObject, eventdata, handles)



function mainAxUpdate(hObject, eventdata, handles)
% hold all
handles.figure1.CurrentAxes = handles.MainPlot;
if ~isempty(handles.xvals) && ~isempty(handles.yvals)
    if strcmp(handles.axislabelX, handles.axislabelY)
        histogram(handles.xvals(:), 50, 'Normalization', 'probability')
    else
%         hashtml = ~cellfun(@isempty, handles.allparams.html);
%         hashtml = hashtml(:);
        if ~isempty(strfind(handles.axislabelY, 'FC'))
            basecolor = [.5 .5 .5];
            basesize = 10;
        else
            basecolor = [1 1 1];
            basesize = 30;
        end
        
        hasmut = handles.allparams.mutmat==1;
        hasmut = hasmut(:);
        
        sizes = ones(length(hasmut), 1).*basesize;
        sizes(hasmut==1, :) = repmat(100, sum(hasmut==1), 1);
        
%         color = ones(length(hasmut), 3);
        color = repmat(basecolor, length(hasmut), 1);
        color(hasmut==1, :) = repmat(handles.colors.bluel, sum(hasmut==1), 1);
        
        scatter(handles.xvals(:), handles.yvals(:), sizes, color);
%         x = handles.xvals(hashtml);
%         y = handles.yvals(hashtml);
%         scatter(x,y);
    end
    
    set(gca, 'Color', [0,0,0], 'XColor',[1 1 1], 'YColor',[1 1 1],  'box','off', 'HitTest','off')
    %handles.MainPlot.Children.ButtonDownFcn = @(hObject, eventdata, handles) MainPlot_ButtonDownFcn;
    guidata(hObject, handles);
%     handles.MainPlot.Children
    for i = 1:length(handles.MainPlot.Children)
        handles.MainPlot.Children(i).ButtonDownFcn = {@MainPlot_ButtonDownFcn,handles};
    end
    
    ylabel(handles.axislabelY, 'color', [1 1 1]);
    xlabel(handles.axislabelX, 'color', [1 1 1]);
%     datacursormode on
%     dcm = datacursormode(handles.figure1);
%     set(dcm,'UpdateFcn',{@MainPlot_ButtonDownFcn,hObject, eventdata, handles})
end

   
% --- Executes on mouse press over axes background.
function MainPlot_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to MainPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

htmlcell = handles.allparams.html;

cP = handles.MainPlot.CurrentPoint;
cx = cP(1,1);
cy = cP(1,2);

x = handles.xvals(:);
y = handles.yvals(:);

diffValues = ((cx-x)./max(x)).^2 + ((cy-y)./max(y)).^2; % note the normalization is required, esp since x and y can have drastically different magnitudes
[~,minValue] = min(diffValues);
[rowClick, colClick] = ind2sub(size(handles.xvals), minValue);
strainclick = handles.allparams.strains{colClick};
readsclick = handles.allparams.summaryReads(colClick);
alignclick = handles.allparams.summaryAlign(colClick);

clickvalues = num2cell([cx, cy, minValue, colClick, rowClick]);
clickvalues = cellfun(@num2str, clickvalues, 'UniformOutput', 0)';

clickvalues = [clickvalues; num2str(readsclick, '%10.3e'); num2str(alignclick)];

handles.infoTable.Data = clickvalues;
handles.infoTable.RowName = {'x value', 'y value', 'index', 'column', 'row', 'reads', 'align%'};
handles.strainName.String = strainclick;

url = htmlcell{minValue};
if ~isempty(url)
    web(url)
end



function param = getvalue(hObject, eventdata, handles, val)
selection = handles.listboxparams{val};
disp(selection)
% if ~strcmp(selection, 'mutation') && ~strcmp(selection, 'strain')
    param = handles.allparams.(['param_', selection]);
%     param = zeros(size(handles.allparams.mutmat));
%     for i = 1:size(handles.allparams.mutmat,1)
%         for j = 1:size(handles.allparams.mutmat,2)
%             entry = handles.allparams.(['param_', selection]);
%             entry = entry(i,j);
%             if isempty(entry)
%                 entry = NaN;
%             end
%             param(i,j) = entry;
%         end
%     end
% end


% --- Executes on button press in clearbutton.
function clearbutton_Callback(hObject, eventdata, handles)
handles.figure1.CurrentAxes = handles.MainPlot;
cla;
% hObject    handle to clearbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over buttonHoldAll.
function buttonHoldAll_Callback(hObject, eventdata, handles)
handles.figure1.CurrentAxes = handles.MainPlot;
hold all

% --- Executes on button press in buttonholdoff.
function buttonholdoff_Callback(hObject, eventdata, handles)
handles.figure1.CurrentAxes = handles.MainPlot;
hold off


% --------------------------------------------------------------------
function infoTable_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to infoTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
