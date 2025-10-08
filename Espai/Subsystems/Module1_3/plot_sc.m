function [] = plot_sc(Bodies)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
figure; hold on; grid on; axis equal;
xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
title('Distribución de subcuerpos y CG global');
colors = lines(length(Bodies));

for i = 1:length(Bodies)
    plot3(Bodies(i).r_cg(1), Bodies(i).r_cg(2), Bodies(i).r_cg(3), 'o', ...
          'MarkerSize', 8, 'MarkerFaceColor', colors(i,:), ...
          'DisplayName', Bodies(i).name);
end

plot3(r_cg_total(1), r_cg_total(2), r_cg_total(3), 'kp', ...
      'MarkerSize', 12, 'MarkerFaceColor', 'y', ...
      'DisplayName', 'CG global');

legend show;
end