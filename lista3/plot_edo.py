"""
File: plot_edo.py

Description:
  Plot the data from the Fortran program which solves the EDO for each i iteration.

Dependencies:
- Numpy
- Matplotlib
- OS

Since:
  - 03/2025

Authors:
  - Pedro C. Delbem <pedrodelbem@usp.br>
"""

import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import os

def plot_edo(filename="results.txt", output_folder="images", outputfilename="output.png"):
    """
    Plot the data from the Fortran program which solves the EDO for each i iteration.

    Parameters
    ----------
    filename : str
        The name of the file containing the data to be plotted. Default is "results.txt".
    output_folder : str
        The folder where the image will be saved. Default is "images".
    outputfilename : str
        The name of the file to be saved. Default is "output.png".

    Notes
    -----
    The first line of the file is skipped, assuming it is a header.
    The file must have at least two columns, where the first column is the iteration number and the second column is the value of y at that iteration.
    If there are more than two columns, the third column is assumed to be the second derivative of y.
    """

    #Create the output folder if it does not exist
    os.makedirs(output_folder, exist_ok=True)

    #Read the data from the file, skipping the first line (header)
    data = np.loadtxt(filename, skiprows=1)

    #Detect number of columns
    num_cols = data.shape[1]

    i = data[:, 0]
    y = data[:, 1]

    #Plot
    plt.figure(figsize=(10, 6))
    plt.plot(i, y, label="y(i)", color='blue', marker='o')

    #If there are more than 2 columns, plot the second derivative
    if num_cols >= 3:
        v = data[:, 2]
        plt.plot(i, v, label="y'(i)", color='red', marker='x')

    plt.xlabel("i")
    plt.ylabel("Values")
    plt.title("EDO solution: y(i)" + (" and y'(i)" if num_cols >= 3 else ""))
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    #Save the plot in the specified folder
    file_name = os.path.join(output_folder, outputfilename)
    plt.savefig(file_name)

    #Close the plot
    plt.close()

def main():
    """
    Main function to execute the plotting of EDO data.

    This function prompts the user to input the desired output file name
    for the plot and calls `plot_edo` to read data from 'results.txt',
    generate the plot, and save it in the 'images' folder with the specified
    file name.

    Returns
    -------
    None
    """

    #Prompt the user for the output file name
    print("Enter the output file:")
    outputfilename = input()

    #Plot the EDO data
    plot_edo("results.txt", "images", outputfilename=outputfilename)

if __name__ == '__main__':
    main()
