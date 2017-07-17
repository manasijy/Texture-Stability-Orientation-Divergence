function [rot] = spin2AxisAngle(spin,CS)

 angle= sqrt(spin(1,2)^2 + spin(1,3)^2 + spin(2,3)^2);
 rot.angle = angle/degree;
ax = [spin(2,3)/angle, spin(3,1)/angle, spin(1,2)/angle];
rot.axis = round(Miller(ax(1),ax(2),ax(3) ,CS));

% rot_vector = [rotaxis rotangle];
% rotmatrix = vrrotvec2mat(r);