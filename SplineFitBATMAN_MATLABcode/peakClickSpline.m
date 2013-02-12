function [ppmx, x,y]=peakClickSpline(Xn,ppm,offset, flagmin, flagRatio)
% written 120213 Dr. Jie Hao, Imperial College London

[n,m]=size(Xn);

%click on peaks
[x,y,button] = ginput;
hold on
plot(x,y,'ro')

clicks=length(x);
for i = 1:clicks
    xy(:,i) = [x(i);y(i)];
end
t = 1:clicks;
ts = 1:0.1: clicks;
xys = spline(t,xy,ts);

plot(xys(1,:),xys(2,:),'b-');
if flagRatio
    ppmx = intersectSpecSplineRatio(Xn, offset, x, y, ppm);
else
    ppmx = intersectSpecSpline(Xn, offset, x, y, ppm, flagmin, 0.0005, 0.001);
end

end