clear;
cs = crystalSymmetry('cubic');
ss = specimenSymmetry('mmm');
%% Creating the grid about an ideal orientation For example Brass : 35,45,0.
%%
phi1_id = 35;
phi_id = 45;
phi2_id = 0;
range = 10;
spacing = 2.5;
eul_ini = grid2euler(phi1_id, phi_id, phi2_id,range,spacing);


prompt = 'The VPSC o/p euler angle file name without .txt extension \n';
name = input(prompt);
filename = [name, '.OUT']; 

% The final orientation after lattice spin 
%To be provided from vpsc

eul_fin = import_eul_fin(filename);%, startRow, endRow);



% Calculating the differential

o_ini = orientation('Euler',eul_ini,cs,ss);%.project2FundamentalRegion;
o_fin = orientation('Euler',eul_fin,cs,ss);%.project2FundamentalRegion;

% del_phi1 = (o_fin.phi1 - o_ini.phi1)/degree;
% del_Phi = (o_fin.Phi - o_ini.Phi)/degree;
% del_phi2 = (o_fin.phi2 - o_ini.phi2)/degree;

del_phi1 = (eul_fin(:,1) - eul_ini(:,1))/degree;
del_Phi = (eul_fin(:,2) - eul_ini(:,2))/degree;
del_phi2 = (eul_fin(:,3) - eul_ini(:,3))/degree;
del_eul = [del_phi1, del_Phi, del_phi2]

J=find(del_phi1 > 180); del_phi1(J)=360-del_phi1(J);
J=find(del_Phi > 180); del_Phi(J)=360-del_Phi(J);
J=find(del_phi2 > 180); del_phi2(J)=360-del_phi2(J);
    
q = quiver3(eul_ini(:,1)/degree,eul_ini(:,2)/degree,eul_ini(:,3)/degree,-del_phi1,-del_Phi,-del_phi2);
c = q.Color;
q.Color = 'blue';
q.MaxHeadSize = 0.2;
q.AutoScaleFactor = 2;
q.Marker = 'o';
q.MarkerSize = 3;
%% 
% To plot sections

prompt = 'Which section you want? \n For phi1 section enter 1\n For Phi section enter 2 and for phi2 section enter 3 \n';
section = input(prompt);
secn_set = eul_ini;
secn_set(:,section) = [];
del_eul(:,section) = [];
eul = ["phi1","phi" ,"phi2"];
eul(section) = [];
prompt = 'Please enter the section angle e.g. 45 1\n ';
ang = input(prompt);
J1=find(eul_ini(:,section)== ang); 
if isempty(J1), disp('the section does exist in the list'); 
else
    q1 = quiver(eul_ini(J1,1)/degree,eul_ini(J1,2)/degree,-del_eul(J1,1),-del_eul(J1,2));
end
xlabel(eul(1));
ylabel(eul(2));