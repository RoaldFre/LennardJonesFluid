data = load("/tmp/pressure.txt");

T = 2;

rho = data(:,1);
P = data(:,2);

Pig = T * rho;

plot(rho, P, rho, Pig, 'k');

