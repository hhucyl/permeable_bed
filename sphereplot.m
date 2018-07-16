function sphereplot(px,py,pz,R)
for i=1:numel(px)
    [x,y,z] = sphere(30);
    X = x.*R + px(i);
    Y = y.*R + py(i);
    Z = z.*R + pz(i);
    mesh(X,Y,Z);
    hold on
end
axis equal
end