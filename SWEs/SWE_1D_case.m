% 1d SWE case 
% call function shallow_water_1d and shallow_water_1d_animation
%  Generate the data arrays.
%
clear,clc
  nx = 41;
  nt = 100;
  x_length = 1.0;
  t_length = 0.4;
  g = 9.8;

  [ h_array, uh_array, x, t ] = shallow_water_1d ( nx, nt, x_length, t_length, g );
%%
for i = 1:100
    z = h_array(:,i);
    plot(x,z,'b:','LineWidth',2)
    hold on
    
    title(['t = ' num2str(t(i))])
    xlabel('location (m)')
    ylabel('water level (m)')
    axis([0 1 0 4])
    hold off
    pause(0.1)
end

%%  Create a sequence of JPG images.
%
  shallow_water_1d_animation ( h_array, uh_array, x, t )

  fprintf ( 1, '\n' );
  fprintf ( 1, 'MAKE_ANIMATION:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );