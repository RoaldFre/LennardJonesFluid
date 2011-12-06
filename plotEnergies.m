data = load("/tmp/data.txt");

t = data(:,1);
E = data(:,2);
K = data(:,3);
V = data(:,4);

hold on;
plot(t, E, 'k');
plot(t, K, 'r');
plot(t, V, 'b');
hold off;
