function plot_satellite(shapes, ah)
% plot_satellite plots a composite satellite of various shapes in 3D.
%
% Syntax:
%   plot_satellite(shapes)
%
% Input: 
%   shapes: Array of structures, each containing:
%       - type: 'sphere', 'cylinder', 'box', or 'cube'
%       - dimensions: struct with fields depending on the shape
%       - position: [x, y, z] position (m) of the shape's center of mass
%       - orientation: (Optional) Orientation parameters (if needed)
%   ah: Axes to plot the shapes into. If not given, a new figure is created
%
% This function generates the geometries of each shape, 
% applies rotations and translations, and plots them in a 3D figure. 
% Each shape is colored differently.

% Define a set of colors to use for the shapes
colors = lines(length(shapes));

% Loop over each shape
for i = 1:length(shapes)
    shape = shapes(i);
    
    % Generate the geometry
    switch lower(shape.type)
        case 'sphere'
            [X, Y, Z] = generate_sphere(shape.dimensions.radius);
            geometryType = 'surf';
        case 'cylinder'
            [X, Y, Z] = generate_cylinder( ...
                shape.dimensions.radius, ...
                shape.dimensions.height);
            geometryType = 'surf';
        case 'box'
            [vertices, faces] = generate_box( ...
                shape.dimensions.length, ...
                shape.dimensions.width, ...
                shape.dimensions.height);
            geometryType = 'patch';
        case 'cube'
            [vertices, faces] = generate_cube(shape.dimensions.side);
            geometryType = 'patch';
        otherwise
            error('Unknown shape type: %s', shape.type);
    end
    
    % Apply rotation if any
    if isfield(shape, 'orientation') && ~isempty(shape.orientation)
        
        % Obtain rotation matrix
        if isfield(shape.orientation, 'rotationMatrix')
            R = shape.orientation.rotationMatrix;
        elseif isfield(shape.orientation, 'eulerAngles')
            phi = shape.orientation.eulerAngles(1);
            theta = shape.orientation.eulerAngles(2);
            psi = shape.orientation.eulerAngles(3);
            R = euler_angles_to_rotation_matrix(phi, theta, psi);
        else
            error('Orientation must include a rotation matrix or Euler angles.');
        end

        % Apply rotation
        if strcmp(geometryType, 'surf')
            [X, Y, Z] = rotate_shape(X, Y, Z, R);
        elseif strcmp(geometryType, 'patch')
            vertices = (R * vertices')';
        end

    end
    
    % Apply translation
    if strcmp(geometryType, 'surf')
        X = X + shape.position(1);
        Y = Y + shape.position(2);
        Z = Z + shape.position(3);
    elseif strcmp(geometryType, 'patch')
        vertices(:,1) = vertices(:,1) + shape.position(1);
        vertices(:,2) = vertices(:,2) + shape.position(2);
        vertices(:,3) = vertices(:,3) + shape.position(3);
    end
    
    % Plot the shape
    if strcmp(geometryType, 'surf')
        surf(ah, ...
            X, Y, Z, ...
            'FaceColor', colors(i,:), ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.8);
    elseif strcmp(geometryType, 'patch')
        patch(ah, ...
            'Vertices', vertices, ...
            'Faces', faces, ...
            'FaceColor', colors(i,:), ...
            'EdgeColor', 'none', ...
            'FaceAlpha', 0.8);
    end
end

end

