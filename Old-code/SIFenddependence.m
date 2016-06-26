%tests the sensitivity of our solution to choice of x_end

%choice of values for the endpoint
xend = [50, 100, 200, 400];

%choice of n that will be used
n = [100, 141, 200, 282, 400, 565, 800, 1131, 1600];

KIbend = zeros(length(n),length(xend));
KIIbend = zeros(length(n),length(xend));
KIcomp = zeros(length(n),length(xend));
KIIcomp = zeros(length(n),length(xend));

for i = 1:length(xend)
    for j = 1:length(n)
        [KIbend(j,i), KIIbend(j,i), ~, ~, ~, ~] = nofluidsolution3(n(j), xend(i),0,1,0);
        [KIcomp(j,i), KIIcomp(j,i), ~, ~, ~, ~] = nofluidsolution3(n(j), xend(i),1,0,0);
    end
end