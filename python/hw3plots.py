import numpy as np
from matplotlib import pyplot as plt
import xml.etree.ElementTree as et
from itertools import cycle
import math

# Get data to plot from the XML output file
def get_data(xml_filename):
    output_xml = et.parse(xml_filename).getroot()
    num_cells = int(output_xml.findtext("./data/number_of_cells"))
    cell_length = np.fromstring(output_xml.findtext("./data/cell_length"), sep="\t")
    phi = np.fromstring(output_xml.findtext("./solution/phi"), sep="\t")
    spectral = np.fromstring(output_xml.findtext("./solution/spectral_radius"), sep="\t")
    n_s = len(spectral)
    start = math.floor(1. * n_s / 4.)
    end = math.ceil(3. * n_s / 4.)
    spectral_val = np.mean(spectral[start:end])
    
    x_vals = np.zeros(num_cells * 2)
    phi_ret = np.zeros(num_cells * 2)

    x = 0.0
    for i in range(num_cells):
        phi_ret[0 + 2 * i] = phi[i]
        phi_ret[1 + 2 * i] = phi[i]
        x_vals[0 + 2 * i] = x
        x += cell_length[i]
        x_vals[1 + 2 * i] = x
    
    return x_vals, phi_ret, spectral_val

colors = cycle(['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e'])
shapes = cycle(['^', 'o', 's', '*', '+'])

filenames = ["../problems/p6.2/ss0.{}-out.xml".format(i) for i in [6, 9, 95, 99]]
values = ["0.60", "0.90", "0.95", "0.99"]

plt.figure(0)
for i, filename in enumerate(filenames):
    x, phi, spec = get_data(filename)
    color = next(colors)
    shape = next(shapes)
    label = r"$\Sigma_s = {},\ \rho\approx{}$".format(values[i], round(spec, 3))
    plt.plot(x, phi,
             linestyle='solid', color=color,
             # marker=shape, markerfacecolor='None',
             # markeredgecolor=color, markeredgewidth=1, 
             label=label)
plt.xlabel(r"$x$")
plt.ylabel(r"$\phi{av}$")
plt.grid(True)
plt.legend(fontsize=12, loc='upper left')
plt.tight_layout()
plt.savefig("p6.2.pdf")
plt.close()

values = ["{}.0".format(i) for i in [5,6,7,8,9]]
filenames = ["../problems/p7.1/z{}-out.xml".format(i) for i in [5,6,7,8,9]]

plt.figure(1)
for i, filename in enumerate(filenames):
    x, phi, spec = get_data(filename)
    color = next(colors)
    shape = next(shapes)
    label = r"$Z_0 = {},\ \rho\approx{}$".format(values[i], round(spec, 3))
    plt.plot(x, phi,
             linestyle='solid', color=color,
             # marker=shape, markerfacecolor='None',
             # markeredgecolor=color, markeredgewidth=1, 
             label=label)
plt.xlabel(r"$x$")
plt.ylabel(r"$\phi_{av}$")
plt.grid(True)
plt.legend(fontsize=12, loc='upper left')
plt.tight_layout()
plt.savefig("p7.1.pdf")
plt.close()
    
