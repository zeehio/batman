function [hObject, handles] = deletePlot(hObject, handles)
% written 120213 Dr. Jie Hao, Imperial College London
if (~isempty(handles.h2))
    if( ishandle(handles.h2))
        delete(handles.h2);
        handles.h2 = [];
    end
end

if (~isempty(handles.hspline))
    if( ishandle(handles.hspline))
        delete(handles.hspline);
        handles.hspline = [];
    end
end

if (~isempty(handles.hselect))
    if( ishandle(handles.hselect))
        delete(handles.hselect);
        handles.hselect = [];
    end
end
guidata(hObject, handles);