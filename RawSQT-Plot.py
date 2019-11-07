import sys, getopt
import numpy as np
from matplotlib import rc
from matplotlib import pyplot as plt
from matplotlib import tri as mtri
rc('lines', linewidth=0.75, markersize=3)
rc('font', family='serif', serif='Times', size=10)
rc('text', usetex=True)

wo_index = 0
mode = 'BRDF'
filename = None

try:
    opts_gnu = ['wo-index=', 'mode=']
    opts, args = getopt.getopt(sys.argv[1:], '', opts_gnu)
    for key, val in opts:
        if (key == '--wo-index'):
            wo_index = int(val)
        if (key == '--mode'):
            mode = val
            if (mode != 'BRDF' and 
                mode != 'BTDF'):
                print('mode must be either BRDF or BTDF')
                sys.exit(1)

    if (len(args) != 1):
        print('1 input filename required')
        sys.exit(1)
    else:
        filename = args[0]

except getopt.error as error:
    print(str(error))
    sys.exit(1)

with open(filename, 'rb') as infile:
    lines = infile.readlines()
    wo_thetas = list(map(float, lines[2].decode('utf-8').split(' ')[1:]))
    wo_theta = wo_thetas[wo_index]
    lines = lines[3:]
    lines = [line.decode('utf-8').split(' ') for line in lines]
    lines = [list(map(float, line[1:])) 
        for line in lines if int(line[0]) == wo_index \
        and (float(line[3]) > 0.0) == (mode == 'BRDF')]
    arr = np.array(lines)
    del lines


arr_orig = np.array(arr, copy=True)
try:
    from scipy.spatial import KDTree
    tree = KDTree(arr[:,:2])
    arr_padding = []
    for phi in np.linspace(0, 2 * np.pi, 32):
        x = np.cos(phi)
        y = np.sin(phi)
        fd, fk = tree.query([x, y])
        f = arr[fk, 3]
        arr_padding.append([x, y, 0, f])
finally:
    arr_padding = np.array(arr_padding)
    arr = np.row_stack((arr, arr_padding))

tri = mtri.Triangulation(arr[:,0], arr[:,1])
tri_interp = mtri.LinearTriInterpolator(tri, arr[:,3])
d = np.linspace(-1, 1, 64)
X, Y = np.meshgrid(d, d)

# Full width = 6.875
# Half width = 3.28125 (single column)
figsize = (8, 8)
fig, ax = plt.subplots(
        nrows=1,
        ncols=1,
        sharex=False,
        sharey=False,
        figsize=figsize)
ax.set_title(r'{} for $\theta_o={:.2f}^\circ$'.format(mode, wo_theta * 180 / np.pi))
ax.set_xlim([-1.25, +1.25])
ax.set_ylim([-1.25, +1.25])
ax.set_xticks([-1, 0, +1])
ax.set_yticks([-1, 0, +1])
cf = ax.pcolormesh(X, Y, tri_interp(X, Y), shading='flat', edgecolors='none')
ax.scatter(arr_orig[:,0], arr_orig[:,1], c='k', marker='x')
fig.colorbar(cf, ax=ax)
plt.axis('equal')

fig.tight_layout()
fig.set_size_inches(*figsize)
fig.savefig('Plot.svg')
plt.show()
