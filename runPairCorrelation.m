close all;

numParticles   = '500';
LJcutoff       = '4';
samples        = '400';
samplesLowDens = '100000'; %luckily we have space partitioning...
measWait       = '10';
measInterval   = '0.05';
timestep       = '0.01';
temperature    = '2';
coupling       = '3';


hold on;
lowDensity  = 0.001; %seperately
minDensity  = 0.2;
densityStep = 0.2;
maxDensity  = 1.0;
densities = [minDensity : densityStep : maxDensity];
ndens = length(densities);

lightCol = 0.7; %color of lightest curve


%low density seperately, because it needs more samples
system(['./main ',numParticles,' ',num2str(lowDensity),' -j 200 -m c -l ',LJcutoff,' -s ',samplesLowDens,' -w ',measWait,' -i ',measInterval,' -t ',timestep,' -T ',temperature,' -c ',coupling]);
data = load('/tmp/data.txt');
t = data(:,1);
g = data(:,2);
plot(t, g, 'Color', lightCol * [1, 1, 1], [';',num2str(lowDensity),';']);

%other densities
for i = 1 : ndens
	density = densities(i);
	system(['./main ',numParticles,' ',num2str(density),' -m c -l ',LJcutoff,' -s ',samples,' -w ',measWait,' -i ',measInterval,' -t ',timestep,' -T ',temperature,' -c ',coupling]);
	data = load('/tmp/data.txt');
	t = data(:,1);
	g = data(:,2);
	plot(t, g, 'Color', lightCol * [1,1,1] - lightCol * i / ndens, [';',num2str(density),';']);
	%TODO: get a \rho in there
endfor


filename='pairCorr'
info=['$N$ = ',numParticles,', $\Delta t$ = ',timestep,', $\tau$ = ',coupling,', $T$ = ',temperature,',  $\tau_\mathrm{relax}$ = ',measWait,', $\tau_\mathrm{sample}$ = ',measInterval,', $N_\mathrm{samples}$ = ',samples,' (',samplesLowDens,' for the curve at the lowest density)'];
caption=['Pair correlation function for various densities. Parameters: ',info,'.'];


destdir   = 'latex/images';
relImgDir = 'images';
ylabrule  = '-1.5cm';
xlab      = '$r$ (a.u.)';
ylab      = '$g(r)$';
width     = '1000';
height    = '1000';

axis([0,4,0,1], 'autoy');
makeGraph(filename,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);

hold off;
