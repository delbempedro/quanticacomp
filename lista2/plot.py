"""
File: plot.py

Description:
  Plot the data from the Fortran program for each i interation.

Dependencies:
- Numpy
- Matplotlib
- OS

Since:
  - 03/2025

Authors:
  - Pedro C. Delbem <pedrodelbem@usp.br>
"""

import numpy as np
import matplotlib.pyplot as plt
import os

def plot_from_file(filename, output_folder="images", grafic_file_name="output.txt", grafic_title="", xlabel="x", ylabel="y"):
    """
    Reads a file and creates a simple line plot of the data
    
    Parameters
    ----------
    filename : str
        The name of the file to read
    
    Returns
    -------
    None
    
    Notes
    -----
    The file should have at least two columns, where the first column is the x data and
    the second column is the y data. The file is read using np.loadtxt.
    """

    #Create the output folder if it does not exist
    os.makedirs(output_folder, exist_ok=True)
    
    #Read the file
    data = np.loadtxt(filename)
    
    #Extract the data
    x = data[:, 0]
    y = data[:, 1]
    
    #Create the plot
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, marker='o', linestyle='-', color='b', label='Dados')
    
    #Set the labels
    plt.xlabel('Iterations')
    plt.ylabel('f(x)')
    plt.title(grafic_title)
    plt.grid()
    
    #Save the plot in the specified folder
    file_name = os.path.join(output_folder, grafic_file_name)
    plt.savefig(file_name)

    #Close the plot
    plt.close()

def main():
    """
    Main function to execute the plotting of data from a specified file.

    This function uses `plot_from_file` to read data from 'output.txt'
    and generate a corresponding plot.

    Returns
    -------
    None
    """
    print("Insert data file name:")
    data_file_name = str(input())
    print("Insert grafic file name:")
    grafic_file_name = str(input())
    print("Insert grafic title:")
    grafic_title = str(input())

    plot_from_file(data_file_name, grafic_file_name=grafic_file_name, grafic_title=grafic_title)

if __name__ == '__main__':
    main()