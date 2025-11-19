import os
import glob
import subprocess
import shutil

gstage = 1
glevel = 2
cangle = 5
deltat = 5

ts = 1
nproc = 32
ddims = [512*glevel,1,512*glevel]

dirNameOld = "stage" + str(gstage)
dirNameNew = "trans" + str(gstage)

### step one: gather bck to single processor

# vars
vrs = ['c','press','tpr','u','v','w']

for v in vrs:
# generate input deck
  f = open("input.txt","w+")
  for d in ddims:
    f.write("%i " % d)
  f.write("\n%i\n" % nproc)
  f.write("1\n")
  f.write(v+"\n")
  f.write("0\n")
  f.write("%i " % ts)
  f.write("%i 1\n" % ts)
  f.close()
# run ChangeProc to single processor
  args = ("./ChangeProc")
  subprocess.run(args)

### step two: move old bck files to subfolder and rename new ones

args = ("rm","-r",dirNameOld)
subprocess.run(args)

os.mkdir(dirNameOld)

for file in glob.glob("*.bck"):
  shutil.move(file,dirNameOld)

for file in glob.glob("*.bck2"):
  os.rename(file,file[:-1])

### step three: extend to next gstage

args = ("./Ext","140",str(glevel),str(gstage),str(cangle),str(deltat))
subprocess.run(args)

### step four: scatter bck to org procs

# vars
vrs = ['c','press','tpr','u','v','w']

for v in vrs:
# generate input deck
  f = open("input.txt","w+")
  for d in ddims:
    f.write("%i " % d)
  f.write("1\n")
  f.write("\n%i\n" % nproc)
  f.write(v+"\n")
  f.write("0\n")
  f.write("%i " % ts)
  f.write("%i 1\n" % ts)
  f.close()
# run ChangeProc to single processor
  args = ("./ChangeProc")
  subprocess.run(args)

### step five: move transition bck files to subfolder and rename new ones

args = ("rm","-r",dirNameNew)
subprocess.run(args)

os.mkdir(dirNameNew)

for file in glob.glob("*.bck"):
  shutil.move(file,dirNameNew)

for file in glob.glob("*.bck2"):
  os.rename(file,file[:-1])

