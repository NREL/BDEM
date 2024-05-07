clear all;

% NOTE: Assumes hopper is oriented in the y-direction, with 
%       the hopper length in the x-dir, and depth in the z-dir

% NOTE: Hopper exit is located at y = 0

% Hopper configuration inputs
inc_angle = 45;
bin_height = 0.3;
bin_length = 0.2975;
hopper_depth = 0.0127;
exit_size = 0.0175;
closed_hopper = 1;

% Domain inputs (for centering the hopper)
domain_x = 0.375;
domain_z = 0.015625;

% Calculating hopper points
bin_x_min = (domain_x - bin_length)/2.0;
bin_x_max = domain_x - (domain_x - bin_length)/2.0;
exit_x_min = domain_x/2.0 - exit_size/2.0;
exit_x_max = domain_x/2.0 + exit_size/2.0;
hopper_z_min = (domain_z - hopper_depth)/2.0;
hopper_z_max = domain_z - (domain_z - hopper_depth)/2.0;
hopper_y_min = 0.0;
hopper_y_max = (bin_length/2.0 - exit_size/2.0) / tan(inc_angle*pi/180.0);
bin_y_max = hopper_y_max + bin_height;

fileID = fopen('hopper.stl','w');
fprintf(fileID,'solid Body 1\n');

% Create top of hopper stl elements
fprintf(fileID,'\tfacet normal 0.000000e+00 -1.000000e+00 0.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, bin_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, bin_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, bin_y_max, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal 0.000000e+00 -1.000000e+00 0.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, bin_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, bin_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, bin_y_max, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');



% Create hopper bin front and back stl elements
fprintf(fileID,'\tfacet normal 0.000000e+00 0.000000e+00 -1.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, bin_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, bin_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal 0.000000e+00 0.000000e+00 -1.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, bin_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal 0.000000e+00 0.000000e+00 1.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, bin_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, bin_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal 0.000000e+00 0.000000e+00 1.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, bin_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');



% Create hopper bin left and right side stl elements
fprintf(fileID,'\tfacet normal 1.000000e+00 0.000000e+00 0.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, bin_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, bin_y_max, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal 1.000000e+00 0.000000e+00 0.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, bin_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal -1.000000e+00 0.000000e+00 0.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, bin_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, bin_y_max, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal -1.000000e+00 0.000000e+00 0.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, bin_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');



% Create hopper front and back stl elements
fprintf(fileID,'\tfacet normal 0.000000e+00 0.000000e+00 -1.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_min, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal 0.000000e+00 0.000000e+00 -1.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_min, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal 0.000000e+00 0.000000e+00 -1.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_min, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_min, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal 0.000000e+00 0.000000e+00 -1.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_min, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal 0.000000e+00 0.000000e+00 1.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_min, hopper_z_min);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal 0.000000e+00 0.000000e+00 1.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_min, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal 0.000000e+00 0.000000e+00 1.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_min, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_min, hopper_z_min);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal 0.000000e+00 0.000000e+00 1.000000e+00\n');
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_min, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');



% Create hopper left and right inclined wall stl elements

% First get normal vectors
p1 = [bin_x_min, hopper_y_max, hopper_z_min];
p2 = [exit_x_min, hopper_y_min, hopper_z_min];
p3 = [bin_x_min, hopper_y_max, hopper_z_max];
normal_left = cross(p3-p1, p2-p1);
nleft = normal_left/norm(normal_left);

p4 = [bin_x_max, hopper_y_max, hopper_z_min];
p5 = [exit_x_max, hopper_y_min, hopper_z_min];
p6 = [bin_x_max, hopper_y_max, hopper_z_max];
normal_right = cross(p5-p4, p6-p4);
nright = normal_right/norm(normal_right);


fprintf(fileID,'\tfacet normal %.8f %.8f %.8f\n', nleft(1), nleft(2), nleft(3));
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_min, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal %.8f %.8f %.8f\n', nleft(1), nleft(2), nleft(3));
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_min, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_min, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_min, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal %.8f %.8f %.8f\n', nright(1), nright(2), nright(3));
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, hopper_y_max, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_min, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

fprintf(fileID,'\tfacet normal %.8f %.8f %.8f\n', nright(1), nright(2), nright(3));
fprintf(fileID,'\t\touter loop\n');
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', bin_x_max, hopper_y_max, hopper_z_max);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_min, hopper_z_min);
fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_min, hopper_z_max);
fprintf(fileID,'\t\tendloop\n');
fprintf(fileID,'\tendfacet\n');

if closed_hopper == 1
    % Create hopper exit cover stl elements
    fprintf(fileID,'\tfacet normal 0.000000e+00 1.000000e+00 0.000000e+00\n');
    fprintf(fileID,'\t\touter loop\n');
    fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_min, hopper_z_min);
    fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_min, hopper_z_min);
    fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_min, hopper_z_max);
    fprintf(fileID,'\t\tendloop\n');
    fprintf(fileID,'\tendfacet\n');

    fprintf(fileID,'\tfacet normal 0.000000e+00 1.000000e+00 0.000000e+00\n');
    fprintf(fileID,'\t\touter loop\n');
    fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_min, hopper_y_min, hopper_z_max);
    fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_min, hopper_z_min);
    fprintf(fileID,'\t\t\tvertex %.8f %.8f %.8f\n', exit_x_max, hopper_y_min, hopper_z_max);
    fprintf(fileID,'\t\tendloop\n');
    fprintf(fileID,'\tendfacet\n');
end

fprintf(fileID,'endsolid Body 1');