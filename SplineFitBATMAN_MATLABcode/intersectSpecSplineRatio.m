function ppmx = intersectSpecSplineRatio(Xn, offset, x, y, ppm)
% written 120213 Dr. Jie Hao, Imperial College London
% modified 21/05/2013
[n,m]=size(Xn);
offy=[0:n-1].*offset;
off = repmat(offy',1,m);
D = Xn + off;
[y iy]=sort(y);
x=x(iy);
t = 1:length(x);
xy = [x,y];

ppmx = zeros(2,size(D,1));

for i = 1:length(x)
    [~, xid]= min(abs(ppm-x(i)));
    [~, yid(i,1)] = min(abs(D(:,xid) - y(i)));
    if (i == 1 || i == length(x))
        [~,fl] = min(abs(D([1,end],xid)-y(i)));
        if (fl == 1)
            yid(i,1) = 1;
        else
            yid(i,1) = size(D,1);
        end
    end
end

ts = [];
for i = 1:length(x)-1
    tmp = i:1/(yid(i+1)-yid(i)):i+1;
    ts = [ts, tmp];
end
ts = unique(ts);
ppmx = spline(t,xy',ts);
