# Libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap

# The file name to read in
filename = 'output/control.Exp40.P4.E3.F1.1970-1974.convergence.newclim'

pltvar  = 0 # (0=Tsurf, 1=Tatmos, 2=Tocean, 3=spec. humidity, 4=albedo, 5=precipitation (mm/day), 6=precip (kg m-2 s-1), 7=eva (kg m-2 s-1), 8=crcl (kg m-2 s-1))
plttime = 0 # (0=First January mean, 1=First Feb. mean, ..., 12=second Jan. mean, ..., -1=last Dec. mean)

# Dimensions
nvar  = 9
ntime = 5*12 # 5 years
nx    = 96
ny    = 48

with open(filename, 'rb') as f:
    data = np.fromfile(f, dtype=np.float32)
    array = np.reshape(data, [ntime,nvar,ny,nx])

lats = np.arange(-88.125,88.125+3.75,3.75)
lons = np.arange(1.875,358.125+3.75,3.75)
lons, lats = np.meshgrid(lons, lats)
fig = plt.figure()
m = Basemap(projection='moll',lon_0=180)
x, y = m(lons,lats)
m.contourf(x,y,array[plttime,pltvar,:,:])
plt.title('Tsurf')
m.drawcoastlines()
plt.colorbar()
plt.show()
