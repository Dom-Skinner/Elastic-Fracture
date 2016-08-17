load n465_data
nv = numel(v);

plot(x(1:t-1),hprime_data(nv+1:nv+t-1,1,1), ...
    x(1:t-1),hprime_data(nv+1:nv+t-1,2,1),...
    x(1:t-1),hprime_data(nv+1:nv+t-1,3,1))
xlabel('x')
ylabel('sqrt(x)h''')
legend({['\lambda = ', num2str(lambda(1))], ...
    ['\lambda = ', num2str(lambda(2))],['\lambda = ', num2str(lambda(3))]},...
    'Location','NorthWest')
