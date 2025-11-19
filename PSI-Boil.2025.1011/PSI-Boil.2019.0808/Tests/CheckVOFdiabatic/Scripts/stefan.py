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

A = 6.695e-2 
B = 0.025
C = 0.597*2030

def process(path,flag=True):
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
    data['x'][i] = float(item)

  f.close()
  return data;

def correct(data):
  for i, item in enumerate(data['t']):
    data['t'][i] = 1e3*float(item)

  for i, item in enumerate(data['x']):
    data['x'][i] = 1e6*float(item)

  return data;

def l2(data,mult=1):
  err = 0
  for i, item in enumerate(data['t']):
    anr = 2*A*np.sqrt(B*float(item)/C)
#err = err + 1/mult*((float(data['x'][i]) - anr)/anr)**2
    err = err + 1/mult*((float(data['x'][i]) - anr))**2

  return np.sqrt(err)

vof_l2 = []

vof2 = process('lvl2/xmn.txt')
t0_vof2 = (vof2['x'][0]/2/A)**2*C/B;
vof2['t'] = np.array(vof2['t']);
vof2['t'] = vof2['t'] + t0_vof2;
vof_l2.append(l2(vof2,2))
print(2,l2(vof2,2))
correct(vof2)

vof4 = process('lvl4/xmn.txt')
t0_vof4 = (vof4['x'][0]/2/A)**2*C/B;
vof4['t'] = np.array(vof4['t']);
vof4['t'] = vof4['t'] + t0_vof4;
vof_l2.append(l2(vof4,4))
print(4,l2(vof4,4))
correct(vof4)

vof6 = process('lvl6/xmn.txt')
t0_vof6 = (vof6['x'][0]/2/A)**2*C/B;
vof6['t'] = np.array(vof6['t']);
vof6['t'] = vof6['t'] + t0_vof6;
vof_l2.append(l2(vof6,6))
print(6,l2(vof6,6))
correct(vof6)

cutofft = 0.1*1e3
cutoffx = 2e-4*1e6

grid = [1/2,1/4,1/6]

def scaling(data,idx1,idx2):
  print(data[idx1],data[idx2],grid[idx1],grid[idx2])
  print( np.log(data[idx1]/data[idx2]) / np.log(grid[idx1]/grid[idx2]))

scaling(vof_l2,0,1)
scaling(vof_l2,1,2)

lengt = 500
ant = np.linspace(0,cutofft*1e-3,num=lengt);
anr = 2*A*np.sqrt(B*ant/C);

plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif')
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, ax1 = plt.subplots(figsize=[10,7])

fmt = matplotlib.ticker.StrMethodFormatter("{x}")
ax1.xaxis.set_major_formatter(fmt)
ax1.yaxis.set_major_formatter(fmt)
ax1.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax1.tick_params(axis='x',pad=8)

plt.plot(1000*ant,1e6*anr, label='Theoretical');
plt.plot(vof2['t'],vof2['x'], label='VOF');
plt.plot(vof4['t'],vof4['x'], label='VOF (2x ref.)');
plt.plot(vof6['t'],vof6['x'], label='VOF (3x ref.)');
plt.xlabel('Time [ms]');
plt.ylabel('Interface position [$\mu$m]');
#plt.legend(loc='upper left');
plt.xlim([0, cutofft]);
plt.ylim([0, cutoffx]);
#plt.grid(True);
#ax1.legend(loc=2)

# Create a set of inset Axes: these should fill the bounding box allocated to
# them.
ax2left = 92
ax2right = 100
bott = 184
topp = 193
ax2 = plt.axes([ax2left,bott,ax2right-ax2left,topp])

ax2.xaxis.set_major_formatter(fmt)
ax2.yaxis.set_major_formatter(fmt)
ax2.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

ax2.set_xlim(left=ax2left,right=ax2right)
ax2.set_ylim(bottom=bott,top=topp)
# Manually set the position and relative size of the inset axes within ax1
ip = InsetPosition(ax1, [0.4,0.1,0.55,0.4])
ax2.set_axes_locator(ip)
# Mark the region corresponding to the inset axes on ax1 and draw lines
# in grey linking the two axes.
mark_inset(ax1, ax2, loc1=1, loc2=2, fc="none", ec='0.5')

ax2.plot(1000*ant,1e6*anr, label='Theoretical');
plt.plot(vof2['t'],vof2['x'], label='VOF');
plt.plot(vof4['t'],vof4['x'], label='VOF (2x ref.)');
plt.plot(vof6['t'],vof6['x'], label='VOF (3x ref.)');
ax2.legend(loc=4)
ax2.set_yticks(np.arange(bott,topp,4))
ax2.set_xticks(np.arange(ax2left,ax2right,3))
#ax2.set_xticklabels(ax2.get_xticks(), backgroundcolor='w')

leg = ax1.legend(loc=2,prop={'size':17})
for line in leg.get_lines():
    line.set_linewidth(3.0)
for item in ([ax1.title, ax1.xaxis.label, ax1.yaxis.label] +
             ax1.get_xticklabels() + ax1.get_yticklabels()):
    item.set_fontsize(22)

for item in ([ax2.title, ax2.xaxis.label, ax2.yaxis.label] +
             ax2.get_xticklabels() + ax2.get_yticklabels()):
    item.set_fontsize(12)


#plt.show();
fig.savefig("stefan.pdf",bbox_inches='tight')
fig.savefig("stefan.pgf",bbox_inches='tight')
