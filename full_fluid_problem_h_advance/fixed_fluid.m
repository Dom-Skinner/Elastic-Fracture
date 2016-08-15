load n465_data
KII = zeros(2,length(L));
for i = 1:2
    tm = [1,13];
    hprime = hprime_data(numel(v)+1:end,end,tm(i));

    % The hypothesis is that KII varies with L, and basically does not
    % depend at all on what the fluid part is doing

    

    for j = 1:length(L)

        [v,w] = spacing(x,z,L(j),r);
        gprime = h_to_g(x,z,v,w,t,hprime,lambda(end)); % N.B lambda only enters 
        % as a correction to g' at infinity, so really does not make much
        % difference to the equations.
        fprintf('%d%% completed \n',round(100*j/numel(L)));
        KII(i,j) = gprime(1);

    end
end

hold on
plot(L,KII(1,:),'o-',L,KII(1,:),'o-',L,KII0)
xlabel('L')
ylabel('K_{II}')
legend({'Solved from refence h'' with \lambda = 0.9,L=1',...
    'Solved from refence h'' with \lambda = 0.9,L=3.4', ...
    'Full problem solved for KI=0'})
