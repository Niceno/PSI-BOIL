1. PSI-BOIL
  Compile and run
  Output: tracked_regions.txt
          *.dat
          console output

2. Postprocess for tracked_regions.txt
2.1. Compile post_track.f90
  % ifort post_track.f90 -o post_track.exe  (Intel compiler)
  % gfortran post_track.f90 -o post_track.exe  (gnu compiler)

2.2. Execute post_track.exe
% post_track.exe
  input: tracked_regions.txt
  output: cvol_ave.out     (time averaged cvol for selected bubbles)
          cvol_ave_all.out (time averaged cvol for all the bubbles)
          fort.*           (time history of each bubble)
          bubbles-z.png    (gnuplot)

3. Data format of fort.* 
  Time  ID  cvol x y z
  For example:
  0.00000000E+00   1  1236 -0.44911900E-17 -0.89823900E-19  0.10154100E-01
  0.10127000E-03   1  1240 -0.44767100E-17 -0.89534100E-19  0.10207300E-01
