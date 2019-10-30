import subprocess
import os
import shutil
import glob

nproc = 5

if os.path.exists("results_diverging_fwd.txt"):
  os.remove("results_diverging_fwd.txt")

if os.path.exists("diverging_grid.txt"):
  os.remove("diverging_grid.txt")

if os.path.exists("diverging_cfl.txt"):
  os.remove("diverging_cfl.txt")

os.system("rm -r diverging_*")

glevels = [1,2,4]
cfls = [0.02,0.05,0.1,0.3]

for g in glevels:
    print("start-grid",str(g))
    args = ("mpiexec","--use-hwthread-cpus","-n",str(nproc),"./Boil",str(g),"0.2","0")
#subprocess.check_call(args)
    subprocess.run(args)
    strfolder = "diverging_"+str(g)+"_0.2"
    os.mkdir(strfolder)
    for f in glob.glob("./*.dat"):
      shutil.move(f,strfolder)
#subprocess.run(["mv","*dat",strfolder])
#os.system("mv *dat "+strfolder) 
    
if os.path.exists("results_diverging_fwd.txt"):
  os.rename("results_diverging_fwd.txt","diverging_grid.txt") 

for c in cfls:
    print("start-cfl",str(c))
    args = ("mpiexec","--use-hwthread-cpus","-n",str(nproc),"./Boil","2",str(c),"0")
#subprocess.check_call(args)
    subprocess.run(args)
    strfolder = "diverging_2_"+str(c)
    os.mkdir(strfolder)
    for f in glob.glob("./*.dat"):
      shutil.move(f,strfolder)

if os.path.exists("results_diverging_fwd.txt"):
  os.rename("results_diverging_fwd.txt","diverging_cfl.txt")

