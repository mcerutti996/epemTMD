import numpy as np
import matplotlib.pyplot as plt
import MatplotlibSettings

data = np.loadtxt("epemCrossSection.dat")

f, (ax1, ax2) = plt.subplots(2, 1, sharex = "all", gridspec_kw = dict(width_ratios = [1], height_ratios = [3, 1]))
plt.subplots_adjust(wspace = 0, hspace = 0)

ax1.set(ylabel = r"\textbf{Arbitrary units}")
ax1.set_xlim([0.001, 5])
#ax1.set_ylim([0.000005, 1000])
#ax1.set_xscale("log")
#ax1.set_yscale("log")
ax1.plot(data[:,0], data[:,1], "r-", label = r"\textbf{Opposite sign}")
ax1.plot(data[:,0], data[:,2], "b-", label = r"\textbf{Same sign}")
ax1.legend(fontsize = 18, loc = "lower left")

ax2.set(xlabel = r"\textbf{$qT$ [GeV]}")
ax2.set_ylabel(r"\textbf{Ratio to OS}", fontsize = 16)
ax1.set_xlim([0.001, 5])
ax2.set_ylim([0, 1.5])
#ax2.set_xscale("log")
ax2.plot(data[:,0], data[:,1]/data[:,1], "r--")
ax2.plot(data[:,0], data[:,2]/data[:,1], "b--")

plt.savefig("epemCrossSection.pdf")
plt.close()

