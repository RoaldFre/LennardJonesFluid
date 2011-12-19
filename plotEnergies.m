function plotEnergies(filename, caption)
close all;

data = load("/tmp/data.txt");

t = data(:,1);
E = data(:,2);
K = data(:,3);
V = data(:,4);

destdir   = 'latex/images';
relImgDir = 'images';
ylabrule  = '-1.5cm';
xlab      = 'time (a.u.)';
ylab      = 'Energy (a.u.)';
width     = '1000';
height    = '800';


hold on;
plot(t, E, 'k;Total;');
plot(t, K, 'r;Kinetic;');
plot(t, V, 'b;Potential;');
axis([t(1), t(end), 1.3*min([E;K;V]), 1.4*max([E;K;V])]);

makeGraph(filename,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);
hold off;
