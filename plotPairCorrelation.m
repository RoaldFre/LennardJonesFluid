function plotPairCorrelation(filename, caption)
close all;

data = load("/tmp/data.txt");

t = data(:,1);
g = data(:,2);

destdir   = 'latex/images';
relImgDir = 'images';
ylabrule  = '-1.5cm';
xlab      = '$r$ (a.u.)';
ylab      = '$g(r)$';
width     = '1000';
height    = '800';

plot(t, g);
makeGraph(filename,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);
