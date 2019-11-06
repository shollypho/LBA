


%%
ndt = 30;
%t1 pre-t0 t2 post b1 slope/drift
yt = @(x,b1,t1,t2) b1.*(delta(x,t1).*(x-t1));

lbasig = yt(1:100,0.5,30,0);

ramp = lbasig(ndt:end);

figure;plot(lbasig);hold on;plot(ramp)
%%

lag = 10;%lag is a vector of lags
%for i = 1:numel(lag)
t1=lag;%lag(i)
t2=ndt-t1;

regrm = [zeros(1,t1) ramp,ones(1,t2).*ramp(end)];%regrm(i)
plot(regrm)

[rho,pval] = corr(regrm(:), data(:),'type','Spearman');%rho(i) for each lag


