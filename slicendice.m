x = test2.x;
y = test2.y;
z = test2.z;
p = test2.p;
p = reshape(p,25,25,25);
[m,i] = max(p,[],"all");
[xi,yi,zi] = ind2sub([25,25,25],i);
slice(p, [15], [18],[15])