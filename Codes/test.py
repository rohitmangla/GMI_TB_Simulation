import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

alpha = ['ABC', 'DEF', 'GHI', 'JKL']

data = np.random.random((4,4))

fig = plt.figure()
ax = fig.add_subplot(111)
cax = ax.matshow(data, interpolation='nearest', norm=LogNorm(vmin=np.nanmin(data), vmax=np.nanmax(data)))
fig.colorbar(cax)

ax.set_xticklabels(['']+alpha)
ax.set_yticklabels(['']+alpha)

plt.show()