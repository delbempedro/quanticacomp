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
import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
import os
import re

def extract_data(filename):
    """
    Extract the relevant data from the given file.

    The file is expected to be as produced by the Fortran program,
    with a header containing the method name, the initial guess,
    the pre-initial guess (if any), the function expression and
    the data points.

    Returns a tuple containing the method name, the function expression,
    the initial guess, the pre-initial guess (or None if not present)
    and an array of shape (n, 2) containing the iteration number
    and the absolute value of the function value at that iteration.

    Parameters:
    - filename: str, the name of the file to read.

    Returns:
    - method: str, the method name.
    - fx_expr: str, the function expression.
    - initial_guess: float, the initial guess.
    - pre_initial_guess: float or None, the pre-initial guess (if any).
    - data: ndarray of shape (n, 2) containing the iteration number and
      the absolute value of the function value at that iteration.
    """
    #Read the file and extract the relevant data
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    #Extract the method name, initial guess, pre-initial guess (if any),
    method = lines[0].strip()
    initial_guess = float(lines[1].split(':')[1].strip())
    pre_initial_guess = None
    if 'Pre-initial guess' in lines[2]:
        pre_initial_guess = float(lines[2].split(':')[1].strip())
        fx_expr = lines[3].split('=')[1].strip()
        data_start = 4
    else:
        fx_expr = lines[2].split('=')[1].strip()
        data_start = 3

    #Extract the data points
    data = []
    for line in lines[data_start:]:
        match = re.match(r'^\s*(\d+)\s+([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)', line)
        if match:
            data.append((int(match.group(1)), float(match.group(2))))
    
    #Convert the data to a numpy array
    return method, fx_expr, initial_guess, pre_initial_guess, np.array(data)

def plot_convergence(file1, file2, output_folder="images"):
    """
    Plot the convergence of two numerical methods from data files.

    This function reads convergence data from two files, which are expected
    to be formatted with headers and data points as produced by a Fortran program.
    It plots the convergence of each method on a log scale, comparing the 
    absolute values of the function values at each iteration.

    Parameters:
    - file1: str, the filename containing data for the first method.
    - file2: str, the filename containing data for the second method.
    - output_folder: str, optional, the folder where the plot image will be saved.
                     Defaults to 'images'.

    Raises:
    - ValueError: if the function expressions from the two files do not match.

    The plot is saved to the specified output folder, with a filename generated
    based on the function expression and initial guesses.
    """
    #Create the output folder if it does not exist
    os.makedirs(output_folder, exist_ok=True)

    #Extract data from the files
    m1, f1_expr, g1, p1, d1 = extract_data(file1)
    m2, f2_expr, g2, p2, d2 = extract_data(file2)

    #Check if the function expressions are the same
    if f1_expr != f2_expr:
        raise ValueError("Both files must refer to the same function.")

    #Check if the methods are the same
    x1, y1 = d1[:, 0], np.abs(d1[:, 1])
    x2, y2 = d2[:, 0], np.abs(d2[:, 1])

    #Check if the methods are the same
    plt.figure(figsize=(10, 6))
    plt.yscale("log")

    #Set the title and labels
    label1 = f"{m1} | Initial: {g1}" + (f" | Pre: {p1}" if p1 is not None else "")
    label2 = f"{m2} | Initial: {g2}" + (f" | Pre: {p2}" if p2 is not None else "")
    title = f"Convergence for f(x) = {f1_expr}\n{label1} & {label2}"

    #Plot the data
    plt.plot(x1, y1, 'o-', label=label1)
    plt.plot(x2, y2, 'x--', label=label2)

    #Set the labels and title
    plt.xscale("linear")
    plt.yscale("log")
    plt.xlabel("Iteration")
    plt.ylabel("|f(x)| (log scale)")
    plt.title(title)
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend()

    #Generate filename
    f_expr_sanitized = f1_expr.replace(" ", "").replace("^", "")
    g_str1 = str(g1).replace(".", "_")
    p_str1 = str(p1).replace(".", "_") if p1 is not None else "none"
    g_str2 = str(g2).replace(".", "_")
    p_str2 = str(p2).replace(".", "_") if p2 is not None else "none"

    #Create the filename
    filename = f"convergence-{f_expr_sanitized}-g{g_str1}-p{p_str1}_g{g_str2}-p{p_str2}.png"
    filepath = os.path.join(output_folder, filename)

    #Save the plot
    plt.tight_layout()
    plt.savefig(filepath)
    plt.close()
    print(f"Saved plot: {filepath}")

def main():
    """
    Main function to execute the script.

    Plots the convergence of Newton-Raphson and Secant methods for two functions
    and saves the plots to the output folder.
    """
    plot_convergence("newtonraphson1.txt", "secant1.txt")
    plot_convergence("newtonraphson2.txt", "secant2.txt")

#Main function to execute the script
if __name__ == "__main__":
    main()
