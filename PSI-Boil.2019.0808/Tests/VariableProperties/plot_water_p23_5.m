%===============================================================================
%
%  The MATLAB script to plot physical properties in "water_p23_5.txt"
% 
%===============================================================================

% open the file
fid = fopen('water_p23_5.txt');

% skip first two rows
A=fscanf(fid, '%s',[1,8]);
B=fscanf(fid, '%s',[1,8]);

% read the rest
C=fscanf(fid, '%f %f %f %f %f %f %f %f',[8,inf])

% rewind the file (close and open)
fclose(fid);
fid = fopen('water_p23_5.txt');

% plot what you read, takin labels from file
color = ['b','g','r','c','m','y','k'];
for i=1:8
  L=fscanf(fid, '%s', 1); % read a label
  if i>1 
    figure
    plot( C(1,:), C(i,:), color(i-1), 'LineWidth', 2 )
    xlabel('T');
    ylabel(L)
  end
end

%-------------------------------------------------------------------------------
% '$Id: plot_water_p23_5.m,v 1.2 2011/05/28 19:10:38 niceno Exp $'/
%-------------------------------------------------------------------------------
