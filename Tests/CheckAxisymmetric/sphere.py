import subprocess
import os

if os.path.exists("results_sphere.txt"):
  os.remove("results_sphere.txt")

radmults = [0.1,0.2,0.3,0.4,0.5,0.6,0.7]
glevels = [2,4,8,12,16]

for r in radmults:
  for g in glevels:
    print("start",str(g),str(r))
    args = ("./Boil",str(g),str(r))
    subprocess.check_call(args)
#popen = subprocess.Popen(args, stdout=subprocess.PIPE)
#popen.wait()

if os.path.exists("uvw-phi-c-kappa_000000.dat"):
  os.remove("uvw-phi-c-kappa_000000.dat")
if os.path.exists("uvw-phi-c-kappa_000000.plt"):
  os.remove("uvw-phi-c-kappa_000000.plt")

log = open("results_sphere.txt", "r")
print(log.read())
log.close()
