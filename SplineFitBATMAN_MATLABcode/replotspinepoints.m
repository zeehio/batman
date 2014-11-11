function [hObject, handles] = replotspinepoints(hObject, handles)
% written 120213 Dr. Jie Hao, Imperial College London
x = handles.x;
y = handles.y;

clicks=length(x);
set(0,'CurrentFigure', handles.plotFigure);
hold on
handles.hselect = plot(x,y,'rp','markersize',handles.mksize);


for i = 1:clicks
    xy(:,i) = [x(i);y(i)];
end
t = 1:clicks;
ts = 1:0.1: clicks;
xys = spline(t,xy,ts);

handles.hspline = plot(xys(1,:),xys(2,:),'b-');
Xn = handles.nmrdata(handles.dorder,:);
offset = str2num(get(handles.offset,'String'));
if offset == 0
    warndlg('Please input offset > 0.');
    return;
end
ppmerror = str2num(get(handles.ppmRange, 'String'));


if (~isempty(handles.h2))
    if( ishandle(handles.h2))
        delete(handles.h2);
        handles.h2 = [];
    end
end

if handles.flagpoint == 1
    flagmin = 0;
    ppmx = intersectSpecSpline(Xn, offset, x, y, handles.ppm, flagmin, ppmerror);
elseif handles.flagpoint == 2
    flagmin = 1;
    ppmx = intersectSpecSpline(Xn, offset, x, y, handles.ppm, flagmin, ppmerror);
else
    ppmx = intersectSpecSplineRatio(Xn, offset, x, y, handles.ppm);
end
if (size(ppmx,1) ~= 2)
    ppmx = ppmx';
end

if (size(ppmx,2) == length(handles.dorder))
    hold on;
    h2 = plot(ppmx(1,:),ppmx(2,:),'bs');
else
    warndlg('Number of spectra mismatch number of points selected. Please re-select.');
    [hObject, handles] = deletePlot(hObject, handles);
    return;
end

handles.ppmx = ppmx;
handles.h2 = h2;
guidata(hObject, handles);