function ppmx = intersectSpecSpline(Xn, offset, x, y, ppm, minflag, ppmerror)
% written 120213 Dr. Jie Hao, Imperial College London

[n,m]=size(Xn);
offy=[0:n-1].*offset;
[r,c] = size(Xn);
off = repmat(offy',1,c);
D = Xn + off;
[y iy]=sort(y);
x=x(iy);

% ppmx = zeros(2,size(D,1));
ppmx = intersectSpecSplineRatio(Xn, offset, x, y, ppm);
for i = 1:size(ppmx,2)
    ind = find(ppm>=ppmx(1,i)-ppmerror & ppm<= ppmx(1,i)+ppmerror);
    if (~isempty(ind))
        if minflag
            [ny, nyi] = min(D(i,ind));
        else
            [ny, nyi] = max(D(i,ind));
        end
        
        ppmx(1,i) = ppm(ind(nyi));
        ppmx(2,i) = ny;
    end
end