close all;

numParticles = '300';
density      = '0.8';
LJcutoff     = '-1';
samples      = '200';
measWait     = '10';
measInterval = '0.1';
timestep     = '0.001';
temperature  = '2';


hold on;
minDensity  = 0.2;
densityStep = 0.2;
maxDensity  = 1.0;
for density = minDensity : densityStep : maxDensity
	system(['./main ',numParticles,' ',num2str(density),' -m c -r -l ',LJcutoff,' -s ',samples,' -w ',measWait,' -i ',measInterval,' -t ',timestep,' -T ',temperature']);

	data = load('/tmp/data.txt');
	t = data(:,1);
	g = data(:,2);
	plot(t, g, 'Color', 0.3 + 0.7 * (maxDensity - [density, density, density]) / maxDensity, [';',num2str(density),';']);
	%TODO: get a \rho in there
endfor


filename='pairCorr'
info=['$N$ = ',numParticles,', $\\Delta t$ = ',timestep,', $\\tau$ = $10^3 \\Delta t$, $T$ = ',temperature,',  $\\tau_\\mathrm{relax}$ = ',measWait,', $\\tau_\\mathrm{sample}$ = ',measInterval,', $N_\\mathrm{samples}$ = ',samples];
caption=['Pair correlation function for various densities. Parameters: ',info,'.'];


destdir   = 'latex/images';
relImgDir = 'images';
ylabrule  = '-1.5cm';
xlab      = '$r$ (a.u.)';
ylab      = '$g(r)$';
width     = '800';
height    = '600';

axis([0,3,0,1], 'autoy');
makeGraph(filename,caption,destdir,relImgDir,xlab,ylab,ylabrule,width,height);

hold off;
