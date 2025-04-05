"""
File: plot-with-header.py

Description:
  Plot the data from the Fortran program which has a header for each i interation.

Dependencies:
- Numpy
- Matplotlib
- OS

Since:
  - 03/2025

Authors:
  - Pedro C. Delbem <pedrodelbem@usp.br>
"""

#Import necessary libraries
import numpy as np  # type: ignore
import matplotlib.pyplot as plt  # type: ignore
import os

def main():
    """
    Main function to execute the script.

    This function is responsible for calling all other functions in the correct order.
    It creates the output folder if it does not exist, reads all data from files,
    parses the data, plots all data on a single figure, sets plot properties,
    sets the plot title, generates the file name, saves and closes the plot.

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    # Create output folder if it does not exist
    os.makedirs("images", exist_ok=True)

    # List of files and plot styles
    datasets = [
        ("Newton-Raphson 1", "newtonraphson1.txt", "blue"),
        ("Secant 1", "secant1.txt", "red"),
        ("Newton-Raphson 2", "newtonraphson2.txt", "green"),
        ("Secant 2", "secant2.txt", "orange")
    ]

    #Initialize data storage
    data_dict = {}
    initial_guess = None
    pre_initial_guess = None

    # Read all data
    for label, filename, color in datasets:
        with open(filename, 'r') as file:
            lines = file.readlines()

            # Parse initial and pre-initial guess
            for line in lines:
                if "Initial guess" in line:
                    value = float(line.split(":")[1].strip())
                    if initial_guess is None:
                        initial_guess = value
                if "Pre-initial guess" in line:
                    value = float(line.split(":")[1].strip())
                    if pre_initial_guess is None:
                        pre_initial_guess = value

            # Parse iteration data
            data = []
            for line in lines:
                parts = line.strip().split()
                if len(parts) == 2:
                    try:
                        data.append([int(parts[0]), float(parts[1])])
                    except ValueError:
                        continue
            data = np.array(data)
            if data.size > 0:
                data_dict[label] = (data[:, 0], np.abs(data[:, 1]), color)

    #Plot all data on a single figure
    plt.figure(figsize=(10, 6))
    for label, (x, y, color) in data_dict.items():
        plt.plot(x, y, marker='o', linestyle='-', color=color, label=label)

    #Set plot properties
    plt.yscale("log")
    plt.xlabel("Iterations")
    plt.ylabel("|f(x)|")

    #Set plot title
    title_guess = f"Initial Guess = {initial_guess:.2f}"
    if pre_initial_guess is not None:
        title_guess += f", Pre-initial Guess = {pre_initial_guess:.2f}"
    plt.title(f"Convergence of Newton-Raphson and Secant Methods\n{title_guess}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    # Generate file name
    initial_str = str(initial_guess).replace('.', '')
    preinitial_str = str(pre_initial_guess).replace('.', '') if pre_initial_guess is not None else "none"
    filename = f"images/convergence-initial-{initial_str}-preinitial-{preinitial_str}.png"

    # Save and close plot
    plt.savefig(filename)
    plt.close()

#Main function to execute the script
if __name__ == '__main__':
    main()