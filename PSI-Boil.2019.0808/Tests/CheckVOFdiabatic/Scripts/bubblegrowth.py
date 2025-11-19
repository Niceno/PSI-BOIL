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


totalvol = 5.96046e-17*48**3
totalvol = 7.45058e-18*96**3

def process(path,flag=True):
  f = open(path, 'rt')
  reader = csv.reader(f, delimiter=' ')

  lineData = list()

  cols = ['aux','t','aux2','x']
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

  #convert strings to floats
  for i, item in enumerate(data['t']):
    data['t'][i] = float(item)

  for i, item in enumerate(data['x']):
    data['x'][i] = float(item)

  f.close()
  return data;

def processtot(path,flag=True):
  f = open(path, 'rt')
  reader = csv.reader(f, delimiter=' ')

  lineData = list()

  cols = ['aux','t','x']
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

  #convert strings to floats
  for i, item in enumerate(data['t']):
    data['t'][i] = float(item)

  for i, item in enumerate(data['x']):
    data['x'][i] = (6*(totalvol-float(item))/np.pi)**(1/3)

  f.close()
  return data;

bg = 4.063399777477446
deltat = 1.25;
rhov =0.597; 
rhol =958.4;
L =2258e3;
cpl =4215.9;
cpv =2030;
lambdal = 0.679
kappa = lambdal/cpl/rhol

cutofft = 0.0012
cutoffx = 0.000125

lengt = 500
ant = np.linspace(0,cutofft,num=lengt);
anr = 2*bg*np.sqrt(kappa*ant);

start = 3062
vof2_ext = processtot('lvl2/tot.txt')
t0_vof2_ext = (vof2_ext['x'][start]/2/bg)**2/kappa;
vof2_ext['t'] = np.array(vof2_ext['t']);
vof2_ext['t'] = vof2_ext['t'] + t0_vof2_ext - vof2_ext['t'][start];

start = 8661
vof4_ext = processtot('lvl4/tot.txt')
t0_vof4_ext = (vof4_ext['x'][start]/2/bg)**2/kappa;
vof4_ext['t'] = np.array(vof4_ext['t']);
vof4_ext['t'] = vof4_ext['t'] + t0_vof4_ext - vof4_ext['t'][start];

start = 15911
vof6_ext = processtot('lvl6/tot.txt')
t0_vof6_ext = (vof6_ext['x'][start]/2/bg)**2/kappa;
vof6_ext['t'] = np.array(vof6_ext['t']);
vof6_ext['t'] = vof6_ext['t'] + t0_vof6_ext - vof6_ext['t'][start];

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

plt.plot(ant,anr, label='theor');
plt.plot(vof2_ext['t'],vof2_ext['x'], label='vof2_ext');
plt.plot(vof4_ext['t'],vof4_ext['x'], label='vof4_ext');
plt.plot(vof6_ext['t'],vof6_ext['x'], label='vof6_ext');
plt.xlabel('Time [s]');
plt.ylabel('Radius [m]');
plt.legend(loc='lower right');
plt.xlim([0, cutofft]);
plt.ylim([0, cutoffx]);
plt.grid(True);

leg = ax1.legend(loc=4,prop={'size':17})
for line in leg.get_lines():
    line.set_linewidth(3.0)
for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
    item.set_fontsize(22)


#plt.show();

fig.savefig("bubblegrowth.pdf",bbox_inches='tight')
fig.savefig("bubblegrowth.pgf",bbox_inches='tight')

