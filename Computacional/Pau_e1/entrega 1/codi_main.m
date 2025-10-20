%
clear all
clc

%----input data----
L = 4; %m
H = 2; %m
Tin = 273+20; %K
Pin = 101325; %Pa
R = 287;
Vin = 15; %m/s

N=200;
M=200;
delta=1e-5;
PSIstart = 1;
acabar = 0;

radius = 0.5; %m (cylinder radius)

cylinder_exists = 1;   % if 1 then it exists. if 0 then it don't exist
rotating_cylinder = 1; % if 1 then it rotate. if 0 then it don't rotate

%----previous calculations----
RHOin = Pin/(Tin*R); %kg/m3

deltaX = L/N;
deltaY = H/M;

Xcv = zeros(N+1,1);
Xp = zeros(N+2,1);
Ycv = zeros(M+1,1);
Yp = zeros(M+2,1);

for i=1:1:N+1
    Xcv(i,1) = deltaX*(i-1);
end

for j=1:1:M+1
    Ycv(j,1) = deltaY*(j-1);
end

Xp(1,1) = 0;
Xp(N+2,1) = L;
for i=2:1:N+1
    Xp(i,1) = deltaX/2 + deltaX*(i-2);
end

Yp(1,1) = 0;
Yp(M+2,1) = H;
for j=2:1:M+1
    Yp(j,1) = deltaY/2 + deltaY*(j-2);
end

%----initial map----
RHO = zeros(N+2,M+2);
RHO(:,:) = RHOin;
PSI = zeros(N+2,M+2);
PSIestimate = zeros(N+2,M+2);
PSIestimate(:,:) = PSIstart;

aE = zeros(N+2,M+2);
aW = zeros(N+2,M+2);
aS = zeros(N+2,M+2);
aN = zeros(N+2,M+2);
bp = zeros(N+2,M+2);

%----calculation of the cylinder-----
solid_region = zeros(N+2,M+2);
cylinderCENTRE=[L/2 H/2];
if cylinder_exists==1
    for i=1:1:N+2
        for j=1:1:M+2
            distTOcylinderCENTRE = sqrt( (Xp(i,1)-cylinderCENTRE(1,1))^2 + (Yp(j,1)-cylinderCENTRE(1,2))^2 );
            if distTOcylinderCENTRE <= radius
                solid_region(i,j)=1;
                RHO(i,j)=1e-10;
            end
        end
    end
end

%----evaluation of coefficients----
%bottom
aE(:,1)=0; aW(:,1)=0; aN(:,1)=0; aS(:,1)=0; bp(:,1)=0; aP(:,1)=1;
PSIestimate(:,1)=0;
%top
aE(:,M+2)=0; aW(:,M+2)=0; aN(:,M+2)=0; aS(:,M+2)=0; bp(:,M+2)=Vin*H; aP(:,M+2)=1;
PSIestimate(:,M+2)=bp(:,M+2);
%inlet
aE(1,:)=0; aW(1,:)=0; aN(1,:)=0; aS(1,:)=0; aP(1,:)=1;
for j=1:1:M+2
    bp(1,j) = Vin*Yp(j,1);
    PSIestimate(1,j)=bp(1,j);
end
%internal
for i=2:1:N+1
    for j=2:1:M+1
        dPE(i,j) = deltaX;
        dPe = deltaX/2;
        dEe = deltaX/2;

        dPW(i,j) = deltaX;
        dPw = deltaX/2;
        dWw = deltaX/2;

        dPN(i,j) = deltaY;
        dPn = deltaY/2;
        dNn = deltaY/2;

        dPS(i,j) = deltaY;
        dPs = deltaY/2;
        dSs = deltaY/2;

        if i==2 && j==2
            dPW(i,j) = deltaX/2;
            dPw = deltaX/2;
            dWw = 0;
            dPS(i,j) = deltaY/2;
            dPs = deltaY/2;
            dSs = 0;
        elseif i==N+1 && j==2
            dPS(i,j) = deltaY/2;
            dPs = deltaY/2;
            dSs = 0;
            dPE(i,j) = deltaX/2;
            dPe = deltaX/2;
            dEe = 0;
        elseif i==2 && j==M+1
            dPW(i,j) = deltaX/2;
            dPw = deltaX/2;
            dWw = 0;
            dPN(i,j) = deltaY/2;
            dPn = deltaY/2;
            dNn = 0;
        elseif i==N+1 && j==M+1
            dPN(i,j) = deltaY/2;
            dPn = deltaY/2;
            dNn = 0;
            dPE(i,j) = deltaX/2;
            dPe = deltaX/2;
            dEe = 0;
        elseif i==2 
            dPW(i,j) = deltaX/2;
            dPw = deltaX/2;
            dWw = 0;
        elseif i==N+1
            dPE(i,j) = deltaX/2;
            dPe = deltaX/2;
            dEe = 0;
        elseif j==2
            dPS(i,j) = deltaY/2;
            dPs = deltaY/2;
            dSs = 0;
        elseif j==M+1
            dPN(i,j) = deltaY/2;
            dPn = deltaY/2;
            dNn = 0;
        end

        RHOref_RHOe(i,j) = dPE(i,j)/((dPe/(RHOin/RHO(i,j)))+(dEe/(RHOin/RHO(i+1,j))));
        RHOref_RHOw(i,j) = dPW(i,j)/((dPw/(RHOin/RHO(i,j)))+(dWw/(RHOin/RHO(i-1,j))));
        RHOref_RHOn(i,j) = dPN(i,j)/((dPn/(RHOin/RHO(i,j)))+(dNn/(RHOin/RHO(i,j+1))));
        RHOref_RHOs(i,j) = dPS(i,j)/((dPs/(RHOin/RHO(i,j)))+(dSs/(RHOin/RHO(i,j-1))));


        aE(i,j) = RHOref_RHOe(i,j)*(deltaY/dPE(i,j));
        aW(i,j) = RHOref_RHOw(i,j)*(deltaY/dPW(i,j));
        aN(i,j) = RHOref_RHOn(i,j)*(deltaX/dPN(i,j));
        aS(i,j) = RHOref_RHOs(i,j)*(deltaX/dPS(i,j));
        aP(i,j) = aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j);
        bp(i,j) = 0;
    end
end
%outlet
aE(N+2,:)=0; aW(N+2,:)=1; aN(N+2,:)=0; aS(N+2,:)=0; bp(N+2,:)=0; aP(N+2,:)=1;
PSIestimate(N+2,:)=PSIestimate(N+1,:);
PSI=PSIestimate;

while acabar == 0
    for i=2:1:N+1
        for j=2:1:M+1
            if  cylinder_exists==1 && rotating_cylinder==1 && solid_region(i,j)==1
                PSI(i,j) = (Vin*(H/2))-5;
            elseif cylinder_exists==1 && rotating_cylinder==0 && solid_region(i,j)==1
                PSI(i,j) = Vin*(H/2);
            else
                PSI(i,j)=(aE(i,j)*PSIestimate(i+1,j) + aW(i,j)*PSIestimate(i-1,j) + aN(i,j)*PSIestimate(i,j+1) + aS(i,j)*PSIestimate(i,j-1) + bp(i,j))/aP(i,j);
            end
        end
    end
    PSI(N+2,:)=PSI(N+1,:);
    %PSI(:,8)=Vin*H*0.5;

    if max( abs(PSIestimate-PSI) ) < delta
        acabar = 1;
    else
        PSIestimate=PSI;
    end
end


Vxp = zeros(N+2,M+2);
Vyp = zeros(N+2,M+2);
Vp = zeros(N+2,M+2);


for i=2:1:N+1
    for j=2:1:M+1

        Vxn(i,j) = RHOref_RHOn(i,j)*((PSI(i,j+1)-PSI(i,j))/(dPN(i,j)));
        Vxs(i,j) = -RHOref_RHOs(i,j)*((PSI(i,j-1)-PSI(i,j))/(dPS(i,j)));
        Vye(i,j) = -RHOref_RHOe(i,j)*((PSI(i+1,j)-PSI(i,j))/(dPE(i,j)));
        Vyw(i,j) = RHOref_RHOw(i,j)*((PSI(i-1,j)-PSI(i,j))/(dPW(i,j)));

        Vxp(i,j) = (Vxn(i,j)+Vxs(i,j))/2;
        Vyp(i,j) = (Vye(i,j)+Vyw(i,j))/2;

        Vp(i,j) = sqrt(Vxp(i,j)^2 + Vyp(i,j)^2);

    end
end



for i=1:1:N+2
    for j=1:1:M+2

        if i==1
            Vp(i,j) = Vin;
            Vxp(i,j) = Vin;
            Vyp(i,j) = 0;

        elseif i==N+2
            Vp(i,j) = Vp(i-1,j);
            Vxp(i,j) = Vxp(i-1,j);
            Vyp(i,j) = Vyp(i-1,j);

        elseif j==1
            Vp(i,j) = Vxs(i,j+1);
            Vxp(i,j) = Vxs(i,j+1);
            Vyp(i,j) = 0;

        elseif j==M+2
            Vp(i,j) = Vxn(i,j-1);
            Vxp(i,j) = Vxn(i,j-1);
            Vyp(i,j) = 0;

        elseif solid_region(i,j)==1
            Vp(i,j) = 0;
            Vxp(i,j) = 0;
            Vyp(i,j) = 0;
        end
    end
end



figure(1)
hold on
if cylinder_exists==1
theta = 0:0.01:2*pi;
fill(radius*cos(theta)+cylinderCENTRE(1,1), radius*sin(theta)+cylinderCENTRE(1,2),'w');
end
axis equal
xlim([0 L]);
ylim([0 H]);
Xpos=zeros(N+2,M+2);
for ii=1:1:M+2
    Xpos(:,ii)=Xp(:,1);
end
Ypos=zeros(N+2,M+2);
for jj=1:1:M+2
    Ypos(:,jj)=Yp(jj,1);
end
quiver(Xpos,Ypos,Vxp,Vyp)
title('velocity field');
xlabel('length [m]');
ylabel('heigh [m]');
hold off


Cp = zeros(N+2,M+2);
P = zeros(N+2,M+2);

Cp(:,:)=1-((Vp(:,:).^2)/(Vin^2));
P(:,:) = Cp(:,:)*0.5*RHOin*Vin^2+Pin;
L=0;
D=0;
Circulation = 0;

for i=2:1:N+1
    for j=2:1:M+1

        if solid_region(i+1,j)==1 && solid_region(i-1,j)==0 && solid_region(i,j+1)==0 && solid_region(i,j-1)==0
            D=D+P(i,j)*deltaY;
            Circulation = Circulation + Vye(i,j)*deltaY;
        elseif solid_region(i+1,j)==0 && solid_region(i-1,j)==1 && solid_region(i,j+1)==0 && solid_region(i,j-1)==0
            D=D-P(i,j)*deltaY;
            Circulation = Circulation - Vyw(i,j)*deltaY;
        elseif solid_region(i+1,j)==0 && solid_region(i-1,j)==0 && solid_region(i,j+1)==1 && solid_region(i,j-1)==0
            L=L+P(i,j)*deltaX;
            Circulation = Circulation - Vxn(i,j)*deltaX;
        elseif solid_region(i+1,j)==0 && solid_region(i-1,j)==0 && solid_region(i,j+1)==0 && solid_region(i,j-1)==1
            L=L-P(i,j)*deltaX;
            Circulation = Circulation + Vxs(i,j)*deltaX;
        end

    end
end

Cl=L/(0.5*RHOin*Vin^2*pi*radius^2);
Cd=D/(0.5*RHOin*Vin^2*pi*radius^2);



figure(2)
hold on
contourf(Xpos,Ypos,PSI,10);
colorbar;
if cylinder_exists==1
theta = 0:0.01:2*pi;
fill(radius*cos(theta)+cylinderCENTRE(1,1), radius*sin(theta)+cylinderCENTRE(1,2),'w');
end
axis equal
title('stream lines');
xlabel('length [m]');
ylabel('heigh [m]');
hold off



figure(3)
hold on
contourf(Xpos,Ypos,P,10);
colorbar;
if cylinder_exists==1
theta = 0:0.01:2*pi;
fill(radius*cos(theta)+cylinderCENTRE(1,1), radius*sin(theta)+cylinderCENTRE(1,2),'w');
end
axis equal
title('pressure distribution');
xlabel('length [m]');
ylabel('heigh [m]');
hold off










