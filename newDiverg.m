%% 
% 

clear;
cs = crystalSymmetry('cubic');
ss = specimenSymmetry('mmm');
%% 
% Creating the grid about an ideal orientation
% 
% For example Brass : 35,45,0.
%%
v1= 15:5:55;
v2= 25:5:65;
v3= -20:5:20;
[X, Y,Z] = meshgrid(v1,v2,v3);
%% 
% Now converting the grid to phi1,phi, phi2 coordinate arrays.

phi1 = reshape(X,[],1);
phi  = reshape(Y,[],1);
phi2 = reshape(Z,[],1);
eul_ini = [phi1, phi, phi2];

% To reconvert
% xx1 = reshape(phi1,9,9,9);
% yy2  = reshape(phi,9,9,9);
% zz3 = reshape(phi2,9,9,9);
%% 
% The final orientation after lattice spin 
%%
%To be provided from vpsc
eul_fin = [phi1_f, phi_f, phi2_f];
%% 
% Finding the rotaion vector from euler angles of rotation

% final_eul = [phi1_f phi_f phi2_f];
% rotm = eul2rotm(final_eul,ZXZ);
% omega21_dot = rotm(2,1);
% omega32_dot = rotm(3,2);
% omega13_dot = rotm(1,3);
%% 
% Finding the rotaion vector from initial and final orientations

misorientation = orientation('Euler',eul_fin,cs,ss)/orientation('Euler', eul_ini,cs,ss);
%% 
% Finding the rates of change of euler angles

phi1_dot = omega21_dot -phi2*cosd(phi);
phi_dot = omega32_dot*cosd(phi1) + omega13_dot*sind(phi1);
phi2_dot = (omega32_dot*sind(phi1) - omega13_dot*cosd(phi1))/sind(phi)
%%

p = Euler(new_gmatrix(:,(j+1)));%/degree;
%% 
% delp = p - p0;



o1 = orientation('Euler', p*degree,cs,ss);
o0 = orientation('Euler', p0*degree,cs,ss); 
del_O = Euler(o1)/degree - Euler(o0)/degree;

toc

% oR = fundamentalRegion(o1.CS,o1.SS)
% plot(oR)
% hold on
% % plot the orientation as it is
% plot(o1,'markercolor','b','markerSize',10)
% 
% % plot the orientation within the fundamental zone
% plot(o1.project2FundamentalRegion,'markercolor','r','markerSize',10)
% hold off


% mori = inv(ori1) * ori2

%ori = ori.project2FundamentalRegion(largeGrain.meanOrientation)