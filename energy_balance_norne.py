import rips
import time
import grpc
import math
import os
from operator import itemgetter
from ecl.eclfile import EclFile
from ecl.grid import EclGrid
from ecl.summary import EclSum
import matplotlib.pyplot as plt
from energy_well_norne import energywell

start = time.time()

# Fetch energy balance calculation from energy_well.py
days,energy_balance,energy_external,energy_internal,energy_dissipated = energywell()

folder_name = os.path.dirname(os.path.abspath(__file__))
dirname = os.path.join(folder_name, 'resinsight_props')
if os.path.exists(dirname) is False:
    os.mkdir(dirname)

# Write calculation results in text files for future use
eext_filename = os.path.join(dirname, 'energy_external.txt')
eint_filename = os.path.join(dirname, 'energy_internal.txt')
edis_filename = os.path.join(dirname, 'energy_dissipation.txt')
ebal_filename = os.path.join(dirname, 'energy_balance.txt')
days_filename = os.path.join(dirname, 'days.txt')
with open(eext_filename,'w') as txt_file:
    for ii in range(len(energy_external)):
        txt_file.write('%f\n' % energy_external[ii])

with open(eint_filename,'w') as txt_file:
    for ii in range(len(energy_internal)):
        txt_file.write('%f\n' % energy_internal[ii])

with open(edis_filename,'w') as txt_file:
    for ii in range(len(energy_dissipated)):
        txt_file.write('%f\n' % energy_dissipated[ii])

with open(ebal_filename,'w') as txt_file:
    for ii in range(len(energy_balance)):
        txt_file.write('%f\n' % energy_balance[ii])

with open(days_filename,'w') as txt_file:
    for ii in range(2,len(days)):
        txt_file.write('%f\n' % days[ii])

end = time.time()
print("Time elapsed: ", end - start)
print("Transferred all results back")

# Generate plots
plt.figure(1)
plt.title('Energy Change Rate in Norne Model', fontsize=15)
plt.plot(days,energy_balance,color='k',label='Energy Balance')
plt.plot(days,energy_dissipated,color='r',label='Dissipated Energy',linestyle='--')
plt.plot(days,energy_external,color='g',label='External energy change',linestyle='--')
plt.plot(days,energy_internal,color='b',label='Internal energy change',linestyle='--')
plt.legend(loc='best')
plt.ylabel("Energy change rate (J/s)")
plt.xlabel("Days")

plt.figure(2)
plt.title('External Energy Change Rate in Norne Model', fontsize=15)
plt.plot(days,energy_external,color='g',label='External energy change',linestyle='-')
plt.legend(loc='best')
plt.ylabel("Energy change rate (J/s)")
plt.xlabel("Days")

plt.figure(3)
plt.title('Internal Energy Change Rate in Norne Model', fontsize=15)
plt.plot(days,energy_internal,color='b',label='Internal energy change',linestyle='-')
plt.legend(loc='best')
plt.ylabel("Energy change rate (J/s)")
plt.xlabel("Days")

plt.figure(4)
plt.title('Dissipation Rate in Norne Model', fontsize=15)
plt.plot(days,energy_dissipated,color='r',label='Dissipated energy',linestyle='-')
plt.legend(loc='best')
plt.ylabel("Energy change rate (J/s)")
plt.xlabel("Days")

plt.show()