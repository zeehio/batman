function [offset_out] = stackplot(D,x,lab,s,zlines,ptype)
% [offset] = stackplot(D,x,lab,s,zlines,ptype]) - plot a series of 'spectra', stacked up the page
%
% D = (nXm) of data, each row is a separate spectrum
% x = (1Xm) vector giving x axis values. if == [1] then plots vs. indices.
% lab = (nX1) vector of strings giving labels for spectra.
%	If lab==1, numbers from 1 to n. If lab==0 then no labels.
% s = (nX1) Vector giving y-spacing of plots. 
%	If s=='zero', chooses own spacing 
%	If s not given, chooses own spacing 
% zlines = (1x1) if 1 then plots zero lines (default 0)
% ptype = (string) plot line type specification (optional)
%
% offset = (nX1) vector of calculated offsets for each row of D
%
% written 290999 TMDE
% revised 041199 to include x axis values
% revised 250200 to label spectra if required.
% revised 090102 to not require x argument
% revised 060308 to have separate switch for zero lines. Also removed grid.
% (c) 1999-2008 Dr. Timothy M D Ebbels, Imperial College, London
% modified by 120213 JHao

if (nargin<2) x = 1:size(D,2); end
if (nargin<5) zlines=0; end

tol=0.1;
dylab = 0.5;
[r,c] = size(D);

if (ischar(x)) x = str2num(x); end
if (x==1) x = 1:c; end

if (r==1)
   offset = 0;
else
   if (~exist('s','var') | strcmp(s,'zero'))
      o1 = max(-(diff(D)'));
      o1 = [0,o1 + max(o1)*tol];
      offset=cumsum(o1);
      if(exist('s','var') & s=='zero')
         o2 = max(D');
         o2 = o2(1:r-1);
         o2 = [0,o2 + max(o2)*tol];
         offset = max([o1;o2]);
         offset = cumsum(offset);
      end
   else
      offset = 0:s:(r-1)*s;
   end
end
off = repmat(offset',1,c);
D = D + off;

if (exist('ptype','var'))
   plot(x,D',ptype)
else
   plot(x,D')
end

if zlines
   hold on
   plot(x,off','k:')
   hold off
end

% Labels
%%%%% modified by JHao
axis tight
%%%%%%%%%%

if (exist('lab','var'))
   if iscell(lab) lab = char(lab); end
   if (lab~=0)
      if (lab==1) lab = num2str([1:r]'); end
      a = axis;
      xpos = a(1) + 0.05*(a(2)-a(1));
      xpos = ones(r,1)*xpos;
      ypos = offset + dylab*max(diff(offset));
      text(xpos,ypos,lab);
   end
end

axis tight
% Allocate output

if (nargout>0) offset_out = off; end
