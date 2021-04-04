import subprocess
import os
import shutil
import glob

if os.path.exists("results_circlevortex.txt"):
  os.remove("results_circlevortex.txt")

if os.path.exists("circlevortex_grid.txt"):
  os.remove("circlevortex_grid.txt")

if os.path.exists("circlevortex_cfl.txt"):
  os.remove("circlevortex_cfl.txt")

os.system("rm -r circlevortex_*")

glevels = [1,2,4,8]
cfls = [0.02,0.05,0.1,0.3]

for g in glevels:
    print("start-grid",str(g))
    args = ("mpiexec","--use-hwthread-cpus","-n","20","./Boil",str(g),"0.2","0")
#subprocess.check_call(args)
    subprocess.run(args)
    strfolder = "circlevortex_"+str(g)+"_0.2"
    os.mkdir(strfolder)
    for f in glob.glob("./*.dat"):
      shutil.move(f,strfolder)
#subprocess.run(["mv","*dat",strfolder])
#os.system("mv *dat "+strfolder) 
    
if os.path.exists("results_circlevortex.txt"):
  os.rename("results_circlevortex.txt","circlevortex_grid.txt") 

for c in cfls:
    print("start-cfl",str(c))
    args = ("mpiexec","--use-hwthread-cpus","-n","20","./Boil","4",str(c),"0")
#subprocess.check_call(args)
    subprocess.run(args)
    strfolder = "circlevortex_4_"+str(c)
    os.mkdir(strfolder)
    for f in glob.glob("./*.dat"):
      shutil.move(f,strfolder)

if os.path.exists("results_circlevortex.txt"):
  os.rename("results_circlevortex.txt","circlevortex_cfl.txt")

