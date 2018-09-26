%===============================================================================
%
%  The MATLAB script to plot physical properties in "co2.txt"
% 
%===============================================================================

% open the file
fid = fopen('co2.txt');

% skip first two rows
A=fscanf(fid, '%s',[1,8]);
B=fscanf(fid, '%s',[1,8]);

% read the rest
C=fscanf(fid, '%f %f %f %f %f %f %f %f',[8,inf])

% rewind the file (close and open)
fclose(fid);
fid = fopen('co2.txt');

% plot what you read, takin labels from file
color = ['b','g','r','c','m','y','k'];
for i=1:8
  L=fscanf(fid, '%s', 1); % read a label
  if i>2 
    figure
    plot( C(1,:), C(i,:), color(i-1), 'LineWidth', 2 )
    xlabel('T');
    ylabel(L)
  end
end

%-------------------------------------------------------------------------------
% '$Id: plot_co2.m,v 1.1 2011/05/28 17:55:38 niceno Exp $'/
%-------------------------------------------------------------------------------
