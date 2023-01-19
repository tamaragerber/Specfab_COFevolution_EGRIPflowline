clear all
close all

M = csvread('eigenvalues_downstream.csv',1,0);

x = M(:,1);
y = M(:,2);
lx = M(:,3);
ly = M(:,4);
lz = M(:,5);
dl = M(:,6);
vx = M(:,7);
vy = M(:,8);
exx = M(:,9);
eyy = M(:,10);
ezz = M(:,11);
exy = M(:,12);
eyx = M(:,13);

cumexx = M(:,14);
cumeyy = M(:,15);
cumezz= M(:,16);
cumexy = M(:,17);
cumeyx = M(:,18);

[lat,lon] = psn2ll(x,y);
wgs84 = almanac('earth','wgs84','meters');
dist = distance(lat(1),lon(1),lat,lon,wgs84);
%dist = cumsum(d);

figure()
subplot(2,1,1)
plot(dist./1000,cumsum(exx),'linewidth',2.5,'color',[0.6,0,0.2]);hold on
plot(dist./1000,cumsum(eyy),'linewidth',2.5,'color',[0.6,0.4,0.2]);hold on
plot(dist./1000,cumsum(ezz),'linewidth',2.5,'color',[0.6,0.8,0.4])
plot(dist./1000,cumsum(exy),'linewidth',1,'color',[0.2,0.2,0.2]);hold on
%plot(dist./1000,cumsum(eyx),'linewidth',1,'color',[0.4,0.4,0.4]);hold on
title('cumulated strain')
legend('\epsilon_{xx}','\epsilon_{yy}','\epsilon_{zz}','\epsilon_{xy}')
xlabel('downstream distance [km]')
ylabel('cumulated strain')
set(gca,'fontsize',14)
grid on

subplot(2,1,2)
plot(dist/1000,lx,'linewidth',2.5,'color',[0.6,0,0.2]);hold on
plot(dist/1000,ly,'linewidth',2.5,'color',[0.6,0.4,0.2])
plot(dist/1000,lz,'linewidth',2.5,'color',[0.6,0.8,0.4])
plot(dist/1000,dl,'linewidth',1,'color','k')
legend('\lambda_x','\lambda_y','\lambda_z','\Delta \lambda')
title('COF eigenvalues')
xlabel('downstream distance [km]')
set(gca,'fontsize',14)
grid on

