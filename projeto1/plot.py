#!/usr/bin/env python3
"""
plot.py

Module to plot the convergence results of Numerical methods.

Dependencies:
- Use numpy to read data
- Use matplotlib to plot

Since:
- 04/2025

Authors:
- Pedro C. Delbem. <pedrodelbem@usp.br>
"""
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore

def read_data(file_name):
    data = np.loadtxt(file_name)
    return data[:, 0], data[:, 1]

def main():

    list_of_coord = []
    list_of_label = ["analytic","zero-infinity","infinity-zero"]

    for i in range(len(list_of_label)):
        file_name = "results-"+list_of_label[i]+".txt"
        list_of_coord.append(read_data(file_name))

    print("Insert delta r:")
    delta_r = float(input())

    plt.figure(figsize=(8, 5))

    for i in range(len(list_of_label)):
        x, y = list_of_coord[i]
        label = list_of_label[i]
        plt.plot(x, y, label=label, marker='o', markersize=3, linestyle='-')
    
    plt.xlabel('$r$')
    plt.ylabel('$\hat{\phi}(r)$')
    plt.title('All Solutions')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f'all-plots-{delta_r}.png', dpi=300)

    for i in range(len(list_of_label)):
        x, y = list_of_coord[i]
        label = list_of_label[i]

        plt.figure(figsize=(8, 5))
        plt.plot(x, y, label=label, marker='o', markersize=4, linestyle='-')
        plt.xlabel('$r$')
        plt.ylabel('$\hat{\phi}(r)$')
        plt.title(f'{label}')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(f'plot-{label}-{delta_r}.png', dpi=300)
        plt.close()

if __name__ == "__main__":
    main()
