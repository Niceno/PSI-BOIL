import subprocess
import os

if os.path.exists("results_cone.txt"):
  os.remove("results_cone.txt")

angles = [10,20,30,40,44,46,50,60,70,80]
glevels = [2,4,8,12,16]

for a in angles:
  for g in glevels:
    print("start",str(g),str(a))
    args = ("./Boil",str(g),str(a))
    subprocess.check_call(args)
#popen = subprocess.Popen(args, stdout=subprocess.PIPE)
#popen.wait()

if os.path.exists("uvw-phi-c-kappa_000000.dat"):
  os.remove("uvw-phi-c-kappa_000000.dat")
if os.path.exists("uvw-phi-c-kappa_000000.plt"):
  os.remove("uvw-phi-c-kappa_000000.plt")

log = open("results_cone.txt", "r")
print(log.read())
log.close()
