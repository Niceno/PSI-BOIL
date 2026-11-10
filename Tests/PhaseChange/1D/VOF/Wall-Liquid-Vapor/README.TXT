Stefan problem: Wall-Liquid-Vapor boiling (2026.11.10)

1. Compile
  edit main.cpp: Level = 1, 2,...
2. Run
  $ ./Boil > log.txt
3. Post process
  $ cat log.txt |grep xmin > front.out
  $ gnuplot
    plot sqrt(2.5e-7-6.6917e-8*x),"front.out" u 2:3 w l
4. Sample flow field
    Results/PNG-Level1/*.png

