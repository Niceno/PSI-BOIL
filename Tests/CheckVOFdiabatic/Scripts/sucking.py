#!/usr/bin/python3
from matplotlib import rc
import numpy as np
import matplotlib as mpl
from matplotlib.ticker import FormatStrFormatter
import matplotlib.pyplot as plt
import matplotlib.font_manager
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from mpl_toolkits.axes_grid1.inset_locator import InsetPosition
import csv
import math as m

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit

mpl.use('pgf')

def figsize(scale):
    fig_width_pt = 417.68646                        # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
    "pgf.texsystem": "pdflatex",
    "text.usetex": True,                # use LaTeX to write all text
    "font.family": "sans-serif",
    "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
    "font.sans-serif": ['Computer Modern Sans serif'],
    "font.monospace": [],
    "font.size": 10,
    "axes.labelsize": 10,               # LaTeX default is 10pt font.
    "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
#    "font.size": 20,
#    "axes.labelsize": 20,               # LaTeX default is 10pt font.
#    "legend.fontsize": 16,               # Make the legend/label fonts a little smaller
#    "xtick.labelsize": 16,
#    "ytick.labelsize": 16,
    "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
    "pgf.preamble": [
    r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
    r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
    ]
}
mpl.rcParams.update(pgf_with_latex)


def process(path,flag=True):
  f = open(path, 'rt')
  reader = csv.reader(f, delimiter=' ')

  lineData = list()

  if flag:
    cols = ['x','T']
  else:
    cols = ['x','aux','T']
  #print(cols)

  for col in cols:
    # Create a list in lineData for each column of data.
    lineData.append(list())


  for line in reader:
    for i in range(0, len(lineData)):
      # Copy the data from the line into the correct columns.
      lineData[i].append(line[i])

  data = dict()

  for i in range(0, len(cols)):
    # Create each key in the dict with the data in its column.
    data[cols[i]] = lineData[i]

  #print(data)

  #convert strings to floats
  for i, item in enumerate(data['x']):
    data['x'][i] = float(item)

  for i, item in enumerate(data['T']):
    data['T'][i] = float(item)
  
  f.close()
  return data;

def l2(theor,sim):
  f = interp1d(sim['x'],sim['T'])
  err = 0.0
  for i, item in enumerate(theor['T']):
    im = i-1
    if im == -1:
      im = i
    ip = i+1
    if ip == len(theor['T']):
      ip = i 
    x1 = theor['x'][im]
    x2 = theor['x'][ip]
    err += (item-f(theor['x'][i]))**2*0.5*(x2-x1)
  
  return [f,m.sqrt(err)];

vof = []
theor = process('theor.new')
theor_err = process('theor.err')

vof.append(process('lvl2/temperature.txt',False))
vof.append(process('lvl4/temperature.txt',False))
vof.append(process('lvl6/temperature.txt',False))

grid = [1/2,1/4,1/6]
l2v_vof = []
l2n_vof = []

for i in range(0, len(vof)):
  l2v_vof.append(l2(theor_err,vof[i]))

for item in l2v_vof:
#  print(item[1])
  l2n_vof.append(item[1])

def func(x, a, b):
  return a * x**b 

def scaling(data,idx1,idx2):
  print(data[idx1],data[idx2],grid[idx1],grid[idx2])
  print( np.log(data[idx1]/data[idx2]) / np.log(grid[idx1]/grid[idx2]))
  
scaling(l2n_vof,0,1)
scaling(l2n_vof,1,2)

def correct(data):
  for i, item in enumerate(data['x']):
    data['x'][i] = 1e3*float(item)

  for i, item in enumerate(data['T']):
    data['T'][i] = float(item) - 373.15

  return data;

correct(theor)

correct(vof[0])
correct(vof[1])
correct(vof[2])

plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax1 = plt.subplots(figsize=[10,7])

fmt = matplotlib.ticker.StrMethodFormatter("{x}")
ax1.xaxis.set_major_formatter(fmt)
ax1.yaxis.set_major_formatter(fmt)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.tick_params(axis='x',pad=8)

xminn = 3.3
xmaxn = 4.5
Tminn = -0.1
Tmaxn = 5 

'''
popt_vof,pcov_vof = curve_fit(func,np.array(grid[2:]),np.array(l2n_vof[2:]))
print(popt_vof)
popt_clr,pcov_clr = curve_fit(func,np.array(grid[2:]),np.array(l2n_clr[2:]))
print(popt_clr)

#plt.plot(grid,l2n_vof)
#plt.plot(grid,func(np.array(grid),*popt_vof))
plt.semilogx(grid,l2n_vof)
plt.semilogx(grid,func(np.array(grid),*popt_vof))
plt.semilogx(grid,l2n_clr)
plt.semilogx(grid,func(np.array(grid),*popt_clr))

'''
plt.plot(theor['x'],theor['T'],label='Theoretical')
plt.plot(vof[0]['x'],vof[0]['T'],label='VOF')
plt.plot(vof[1]['x'],vof[1]['T'],label='VOF (2x refined)')
plt.plot(vof[2]['x'],vof[2]['T'],label='VOF (3x refined)')
plt.xlim([xminn,xmaxn]);
plt.ylim([Tminn,Tmaxn]);
plt.xlabel('Distance [mm]');
plt.ylabel('Superheat [K]');

leg = ax1.legend(loc=4,prop={'size':17})
for line in leg.get_lines():
    line.set_linewidth(3.0)
for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
    item.set_fontsize(22)

#plt.show();
fig.savefig("sucking.pdf",bbox_inches='tight')
fig.savefig("sucking.pgf",bbox_inches='tight')
