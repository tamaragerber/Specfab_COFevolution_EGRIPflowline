clear all
close all

fl = shaperead('/home/tamara/Documents/university_copenhagen/PhD_project/Manuscript_NEGIS_fabric/figures/alle_Profile_csv/alle_Profile_csv/D_csv_8/downstream_flowline_rawdata.shp');

X = extractfield(fl,'X');
Y = extractfield(fl,'Y');
vx = extractfield(fl,'vx_clipped');
vy = extractfield(fl,'vy_clipped');
alpha =extractfield(fl,'alpha_FLOW');
dudx =extractfield(fl,'dudx_or');
dudy =extractfield(fl,'dudy_or');
dvdx =extractfield(fl,'dvdx_or');
dvdy =extractfield(fl,'dvdy_or');

x = X;
y = Y;
alpha = -alpha;

[lat,lon] = psn2ll(X,Y);
wgs84 = almanac('earth','wgs84','meters');
dist = distance(lat(1),lon(1),lat,lon,wgs84);

for i = 1:length(lat)
    
    R = [cosd(alpha(i)), -sind(alpha(i));
    sind(alpha(i)), cosd(alpha(i))];


    v_orig = [vx(i);vy(i)];

    v_flow = R*v_orig;

    vx_flow(i) = v_flow(1);
    vy_flow(i) = v_flow(2);

    strain = [dudx(i), 0.5*(dudy(i)+dvdx(i));
            0.5*(dvdx(i)+dudy(i)), dvdy(i)];
        
    spin = [0, 0.5*(dudy(i)-dvdx(i));
            0.5*(dvdx(i)-dudy(i)), 0];
        
        
        strain_flow = R*strain*R';
        
        spin_flow = R*spin*R';
        
     exx_flow(i) = strain_flow(1,1);
     eyy_flow(i) = strain_flow(2,2);
     exy_flow(i) = strain_flow(1,2);
     eyx_flow(i) = strain_flow(2,1);
     
     sxx_flow(i) = spin_flow(1,1);
     syy_flow(i) = spin_flow(2,2);
     sxy_flow(i) = spin_flow(1,2);
     syx_flow(i) = spin_flow(2,1);
      
end

figure()
plot(dist,vx);hold on
plot(dist,vy)
plot(dist,sqrt(vx.^2 + vy.^2),'k');

plot(dist,vx_flow);hold on
plot(dist,vy_flow)
plot(dist, sqrt(vx_flow.^2+vy_flow.^2),'k:');
legend('vx','vy','vmag','vxflow','vyflow','vmag')


figure()
plot(dist,dudx);hold on
plot(dist,dvdy);hold on
plot(dist,dudy);hold on
plot(dist,dvdx)

plot(dist,exx_flow);hold on
plot(dist,eyy_flow);hold on
plot(dist,exy_flow);hold on
plot(dist,eyx_flow);hold on
ezz = 0-exx_flow-eyy_flow;

figure()
subplot(2,1,1)
plot(dist/1000,cumsum(exx_flow));hold on
plot(dist/1000,cumsum(eyy_flow));hold on
plot(dist/1000,cumsum(exy_flow));hold on
plot(dist/1000,cumsum(eyx_flow));hold on
plot(dist/1000,cumsum(ezz))
legend('exx','eyy','exy','eyx','ezz')
xlabel('downstream distance [km]')
ylabel('accumulated strain]')

subplot(2,1,2)
plot(dist/1000,vx_flow);hold on
plot(dist/1000,vy_flow)
plot(dist/1000, sqrt(vx_flow.^2+vy_flow.^2),'k:');
legend('vxflow','vyflow','vmag')
xlabel('downstream distance [km]')
ylabel('velocity [ma^-1]')


v_mag = sqrt(vx_flow.^2+vy_flow.^2);

%% saving

for i = 1:length(fl)-1
    dt(i) = sqrt((x(i+1)-x(i))^2 + (y(i+1)-y(i))^2)/(0.5*(v_mag(i)+v_mag(i+1)));
end
dt(length(fl))=1;

m = [x',y',vx_flow',vy_flow',exx_flow',eyy_flow',ezz',exy_flow',eyx_flow',dt',sxx_flow',syy_flow',sxy_flow',syx_flow']; 
writematrix(m,'strainrates_spin.csv')
