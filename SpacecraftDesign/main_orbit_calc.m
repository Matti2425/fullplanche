%% Init
clearvars
close all

moon.r = 1.737e6;
moon.mu = 4.9028e12;

station.rp = 3.3e6;
station.ra = 7.1e7;

station.a = (station.rp+station.ra)/2;

%% Case 1: Leaving station at apoapsis

hoh1.a = (station.ra + moon.r)/2;
hoh1.ra = station.ra; 

hoh1.dV1 = calc_ellipse_vel(moon.mu, hoh1.a, hoh1.ra) - calc_ellipse_vel(moon.mu, station.a, station.ra);
hoh1.dV2 = 0 - calc_ellipse_vel(moon.mu, hoh1.a, moon.r); % vf = 0;

hoh1.totaldV = hoh1.dV1 + hoh1.dV2;

%% Case 2: Leaving station at periapsis

hoh2.a = (station.rp + moon.r)/2;
hoh2.ra = station.rp; 

hoh2.dV1 = calc_ellipse_vel(moon.mu, hoh2.a, hoh2.ra) - calc_ellipse_vel(moon.mu, station.a, station.rp);
hoh2.dV2 = 0 - calc_ellipse_vel(moon.mu, hoh2.a, moon.r); % vf = 0;

hoh2.totaldV = hoh2.dV1 + hoh2.dV2;
