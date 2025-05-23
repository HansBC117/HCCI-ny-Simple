% Surf_example

% This is just to create an example
% =================================
t1 = linspace(0,2*pi,10);
t2 = linspace(0,2*pi,10);
for i = 1:10
curves(i,1:10) = sin(t2(i))*sin(t1);
end
% =================================

% Grids are made based on vectors. In task 2 the vectors will be
% CAD and injection timing 
[X,Y] = meshgrid(t1,t2);
Z = curves; % In task 2 this will be a matix with HRR

% Interpolation in order to change resolution and get nicer plot	
[Xi,Yi] = meshgrid(linspace(0,2*pi,100),linspace(0,2*pi,100));
Zi = interp2(X,Y,Z,Xi,Yi,'spline');

figure
surf(Xi,Yi,Zi)
caxis([0 1]) % manual control of value range for the color map
view([0 0 1]) % manual change of view angle
shading flat % To remove grid lines 
colorbar 
%xlim([0 2*pi])
%ylim([0 2*pi])

xlabel('Angle 1')
ylabel('Angle 2')
zlabel('Amplitude')

title('Surf example')








% Example call for filtering when 'time' and 'signal' variables are available
% filtered_signal = fftBPfilter(time, signal, [0 4000 1 0], 2000, 'plot_off');





%Hans Code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

[X,Y] = meshgrid(CAD_vector, injection_timings);
Zi = interp2(X,Y,HRR_matrix,Xi,Yi,'spline');

surf(Xi,Yi,Zi)
view([0 0 1])
shading flat
colorbar
xlabel('CAD [°]')
ylabel('Injection timing [CAD BTDC]')
title('Relative HRR map')