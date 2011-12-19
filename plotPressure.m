function plotPressure(filename, caption, T)
close all;

data = load("/tmp/pressure.txt");

rho = data(:,1);
P = data(:,2);

rhoig = [0; 2*max(rho)];
Pig = T * rhoig;

hold on;
plot(rho,   P,   'b;Simulation;');
plot(rhoig, Pig, 'k;Ideal gas;');
legend('location', 'northwest');

destdir   = 'latex/images';
relImgDir = 'images';
ylabrule  = '-1.5cm';
xlab      = '$\rho$ (a.u.)';
ylab      = '$P(\rho)$ (a.u.)';
width     = '1000';
height    = '1000';

axis([0, 1.0*max(rho), 0, 1.0*max([P;Pig])]);

makeGraph(filename,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);
hold off;

