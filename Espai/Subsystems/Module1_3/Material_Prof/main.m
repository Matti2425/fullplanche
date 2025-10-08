%This is my version of the main script to generate the inertia tensor
%Take into account that the positioning of the object is made based on the 
%center of mass

%% GPS Satellite
%Body
shape_body.type = 'box'; 
shape_body.mass = 2000; 
shape_body.dimensions.length = 2;
shape_body.dimensions.width = 2;
shape_body.dimensions.height = 2;
shape_body.position = [0, 0, 0];  
shape_body.orientation = []; 

%Solar Arrays
shape_solarArray1.type = 'box'; 
shape_solarArray1.mass = 50; 
shape_solarArray1.dimensions.length = 4; 
shape_solarArray1.dimensions.width = 0.2;
shape_solarArray1.dimensions.height = 2;
shape_solarArray1.position = [3, 0, 0]; 
shape_solarArray1.orientation = [];


% Solar Array 2
shape_solarArray2.type = 'box'; 
shape_solarArray2.mass = 50; 
shape_solarArray2.dimensions.length = 4; 
shape_solarArray2.dimensions.width = 0.2;
shape_solarArray2.dimensions.height = 2; 
shape_solarArray2.position = [-3, 0, 0]; 
shape_solarArray2.orientation = [];

shapes = [ shape_body;
    shape_solarArray1;
    shape_solarArray2;
    ];

% Plot the satellite
%fh = figure;
hold on;
axis equal;
grid on;
view(3);
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
title('Satellite Geometry');
plot_satellite(shapes, gca);