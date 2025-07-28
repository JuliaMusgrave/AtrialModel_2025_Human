''' 
Short script to help prepare CellML code with piecewise Calcium as an input.
Reads the csv files containing Ca data and outputs a text file suitable for copying into CellML text editor
Useful if data can not be easily captured with a fitted function.
'''
import csv
import numpy as np

ca_time = []
nd_ca = []
d_ca = []
with open('nd_ca.csv', 'r') as file:
    csv_reader = csv.reader(file)
    next(csv_reader)
    for row in csv_reader:
        ca_time.append(float(row[0]))
        nd_ca.append(float(row[1]))

with open('d_ca.csv', 'r') as file:
    csv_reader = csv.reader(file)
    next(csv_reader)
    for row in csv_reader:
        d_ca.append(float(row[1]))

dt = 0.001 # sampling frequency
T = 1 # period
with open('nd_step.txt', 'w') as fp:
    fp.write("Ca = sel\n")
    for t in np.linspace(0, T, int(T/dt)+1):
        ca_i = np.interp(t, ca_time, nd_ca)
        fp.write(f"case t_twitch < {t+dt/2:.4f}{{second}}:\n"
                 f"    {ca_i:.4f}{{micromolar}};\n")
    fp.write("endsel;\n")

with open('d_step.txt', 'w') as fp:
    fp.write("Ca = sel\n")
    for t in np.linspace(0, T, int(T/dt)+1):
        ca_i = np.interp(t, ca_time, d_ca)
        fp.write(f"case t_twitch < {t+dt/2:.4f}{{second}}:\n"
                 f"    {ca_i:.4f}{{micromolar}};\n")
    fp.write("endsel;\n")