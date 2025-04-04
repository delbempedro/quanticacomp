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

import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
import os

def main(filename="results.txt", output_folder="images"):
    """
    Main function to execute the plotting of data from a specified file.

    This function uses `np.loadtxt` to read data from the specified file
    and generate a corresponding plot.

    Parameters
    ----------
    filename : str
        The name of the file to read. Default is "results.txt".
    output_folder : str
        The folder where the plot image will be saved. Default is "images".

    Returns
    -------
    None
    """

    #Create the output folder if it does not exist
    os.makedirs(output_folder, exist_ok=True)

    #ead the data from the file, skipping the first line (header)
    data = np.loadtxt(filename, skiprows=1)

    # Colunas: i, y(i), v(i)
    i = data[:, 0]
    y = data[:, 1]
    v = data[:, 2]

    # Plot
    plt.figure(figsize=(10, 6))
    plt.plot(i, y, label="y(i)", color='blue', marker='o')
    plt.plot(i, v, label="y'(i)", color='red', marker='x')
    plt.xlabel("i")
    plt.ylabel("Values")
    plt.title("EDO solution: y(i) e y'(i)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    #Save the plot in the specified folder
    file_name = os.path.join("images", "results-ex-2.png")
    plt.savefig(file_name)

    #Close the plot
    plt.close()

if __name__ == '__main__':
    main()