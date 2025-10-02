%% 1
%a
axis = [0,0,1];
angle= deg2rad(45);
disp('--------------------------------1:')
disp('--->Rotation matrix: ')
M1 = axisAngle2RotMat(axis,angle);
disp(M1)
disp('Determinante:')
disp(det(M1))
disp('--->Axis-angle:')
disp(axis)
disp(angle)
disp('--->Quaternion:')
q1 = axisAngle2Quat(axis,angle);
disp(q1)
disp('Norm: ')
disp(norm(q1))

%% 2
axis = [1,1,1];
angle= deg2rad(120);
disp('--------------------------------2:')
disp('--->Rotation matrix: ')
M2 = axisAngle2RotMat(axis,angle);
disp(M2)
disp('Determinante:')
disp(det(M2))
disp('--->Axis-angle:')
disp(axis)
disp(angle)
disp('--->Quaternion:')
disp(axisAngle2Quat(axis,angle))
disp(norm(axisAngle2Quat(axis,angle)))

%% 3 

%% 4
axis4 = [1,0,0];
angle4 = deg2rad(180);
disp('--------------------------------4:')
disp('--->Rotation matrix: ')
M4 = axisAngle2RotMat(axis4,angle4);
disp(M4)
disp('Determinante:')
disp(det(M4))
disp('--->Axis-angle:')
disp(axis4)
disp(angle4)
disp('--->Quaternion:')
q4 = axisAngle2Quat(axis4,angle4);
disp(q4)
disp('Norm: ')
disp(norm(q4))


%% 5
disp('--------------------------------5:')
axis5 = [1,0,0];
angle5 = 60;
q5 = axisAngle2Quat(axis5,angle5);
disp('--->Quaternion: ')
disp(q5)
disp('--->Modulus: ')
disp(norm(q5))
