
dXdt = @(y,z) y*z;
dYdt = @(x,z) -2*x*z;
dZdt = @(x,y) x*y;

t0 = [0.5 0 0.5];

figure
hold on
vNow = t0;
vNext = vNow;
lastT = 0;
tStep = 0.0001;
for t=tStep:tStep:1
scatter3(vNow(1),vNow(2),vNow(3));
vNext(1) = vNow(1) + tStep*dXdt(vNow(2),vNow(3));
vNext(2) = vNow(2) + tStep*dYdt(vNow(1),vNow(3));
vNext(3) = vNow(3) + tStep*dZdt(vNow(1),vNow(2));
line([vNow(1);vNext(1)],[vNow(2);vNext(2)],[vNow(3);vNext(3)]);

vNow = vNext;
%pause(1);
drawnow
end

%%
F = @(t,y) [y(2)*y(3); -2*y(1)*y(3); y(1)*y(2)];
y0 = [0 1 0]' + .2*randn(3,1);
y0 = y0/norm(y0);
ode45(F,[0 10],y0)

%%
clearvars
[X,Y] = meshgrid(-2:.1:2);
Z = exp(-((X.^2 + Y.^2)-1).^2);
[DX,DY] = gradient(Z,.1,.1);

figure
subplot(2,2,1)
imagesc(Z)
subplot(2,2,2)
contour(X,Y,Z)
hold on
quiver(X,Y,DX,DY)
hold off
