a = 0.2;
x = -1:0.001:1;
plot(x,((1+x).^a - 1)./x)

fun = @(x) ((1+x).^a - 1)./x;

integral(fun, -1,1)