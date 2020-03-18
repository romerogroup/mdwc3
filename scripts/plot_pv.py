#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import re
import sys

#number of abinit iterations
iterations=int(sys.argv[1])

with open('pressure_volume.mdout', 'r') as f:
    f.readline()  
    lines = f.readlines()
    time = [float(line.split()[3]) for line in lines]

    pressure =  [float(line.split()[4]) for line in lines]
    volume = [float(line.split()[5]) for line in lines]


plt.figure(1)

plt.subplot(211)
plt.plot(time,pressure)
plt.title('Pressure vs Time')
plt.xlabel('Time (fs)')
plt.ylabel('Pressure (hartree/Bohr^3)')

plt.subplot(212)
plt.plot(time,volume)
plt.title('Volume vs Time')
plt.xlabel('Time (fs)')
plt.ylabel('Volume (Bohr^3)')

plt.show()

f.close()