clear
load n995x846-late
tm = [1,2,8];
KII = zeros(numel(tm),length(L));
for i = 1:numel(tm)
    
    hprime = h0_prime(:,tm(i)); % hprime_data(numel(v)+1:end,end,tm(i));

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

plot(l0,KII(1,:),'o-',l0,KII(2,:),'o-',l0,KII(3,:),'o-',l0,KII0)
xlabel('\lambda')
ylabel('K_{II}')
legend({'Solved from refence h'' with L = 0.8',...
    'Solved from refence h'' with L=1', ...
    'Solved from refence h'' with L=2.2', ...
    'Full problem solved for KI=0'})


fid = fopen('fixed-fluid.csv','w');
fprintf(fid,'lambda,     L08,    L10,  L22,   KII0\n');
for j = 1:numel(L)
    fprintf(fid, '%.5e,    %.5e,   %.5e,     %.5e,     %.5e\n', ...
        l0(j), 3*sqrt(2*pi)*KII(1,j), 3*sqrt(2*pi)*KII(2,j), ....
        3*sqrt(2*pi)*KII(3,j), 3*sqrt(2*pi)*KII0(j));
end
fclose(fid);
