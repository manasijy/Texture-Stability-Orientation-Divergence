function [eul_ini] = grid2euler(phi1_id, phi_id, phi2_id,range,spacing)

% grid2euler(phi1_id, phi_id, phi2_id)
% This function creates a 5 degree grid around an ideal orientation. Its
% range is -20 to +20 degrees from the ideal orientation in the fundamental
% euler space. It creates a 729 row vector of the three euler angles i.e.
% a 729X3 matrix. This matrix can be used as an input to simulation tools.

%% Creating the grid about an ideal orientation For example Brass : 35,45,0.


v1= (phi1_id - range):spacing:(phi1_id + range);
v2= (phi_id - range):spacing:(phi_id + range);
v3= (phi2_id - range):spacing:(phi2_id + range);
[X, Y,Z] = meshgrid(v1,v2,v3);

% Now converting the grid to phi1,phi, phi2 coordinate arrays.

phi1 = reshape(X,[],1);
phi  = reshape(Y,[],1);
phi2 = reshape(Z,[],1);
nl = length(phi1);
onem = ones(nl,1);
eul_ini = [phi1, phi, phi2];
mat = [eul_ini onem];
fid= fopen('INPUT.TEX','w');%,'n','US-ASCII');
fprintf(fid,'\n \n');
fprintf(fid,'\nB\t%s\t0',num2str(nl));
dlmwrite('INPUT.TEX',mat,'-append','delimiter',' ','roffset',1);
fclose(fid);
clear mat;
eul_ini = [phi1, phi, phi2]*degree;
