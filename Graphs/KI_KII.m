load n350x873-2t

fid1 = fopen('KI-KII-1.csv','w');
fprintf(fid1,'lambda,  L,     KI,  \n');

fid2 = fopen('KI-KII-2.csv','w');
fprintf(fid2,'lambda,  L,     KII,  \n');

for j = 1:5
subplot(2,2,1)
plot(L, KI(4*j-3,:),'+-');
xlabel('L');
ylabel('K_I');
hold on

for l = 1:numel(L)
    fprintf(fid1, '%.5e,    %.5e,    %.5e \n',...
        lambda(4*j-3), L(l), 3*sqrt(2*pi)*KI(4*j-3,l));
end


subplot(2,2,2)
plot(L, KII(4*j-3,:),'o-');
xlabel('L');
ylabel('K_{II}');
hold on

for l = 1:numel(L)
    fprintf(fid2, '%.5e,    %.5e,    %.5e \n',...
        lambda(4*j-3), L(l), 3*sqrt(2*pi)*KII(4*j-3,l));
end

end
legend({['\lambda=',num2str(lambda(1))],['\lambda=',num2str(lambda(5))],['\lambda=',num2str(lambda(9))],...
    ['\lambda=',num2str(lambda(13))],['\lambda=',num2str(lambda(17))]})
fclose(fid1);
fclose(fid2);




fid1 = fopen('KI-KII-3.csv','w');
fprintf(fid1,'lambda,  L,     KI,  \n');

fid2 = fopen('KI-KII-4.csv','w');
fprintf(fid2,'lambda,  L,     KII,  \n');

for j = 1:4
subplot(2,2,3)
plot(lambda, KI(:,2*j-1),'+-');
xlabel('\lambda');
ylabel('K_I');
hold on

for l = 1:numel(lambda)
    fprintf(fid1, '%.5e,    %.5e,    %.5e \n',...
        lambda(l), L(2*j-1), 3*sqrt(2*pi)*KI(l,2*j-1));
end


subplot(2,2,4)
plot(lambda, KII(:,2*j-1),'o-');
xlabel('\lambda');
ylabel('K_{II}');
hold on

for l = 1:numel(lambda)
    fprintf(fid2, '%.5e,    %.5e,    %.5e \n',...
        lambda(l), L(2*j-1), 3*sqrt(2*pi)*KII(l,2*j-1));
end

end
legend({['L=',num2str(L(1))],['L=',num2str(L(3))],['L=',num2str(L(5))],['L=',num2str(L(7))]})
fclose(fid1);
fclose(fid2);