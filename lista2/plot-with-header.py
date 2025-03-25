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

import numpy as np # type: ignore
import matplotlib.pyplot as plt # type: ignore
import os

def parse_coefficients(coefficients_str):
    """
    Parses the coefficient string and formats it as a LaTeX-formatted polynomial equation.
    """
    coeffs = list(map(float, coefficients_str.split()))
    terms = []
    degree = len(coeffs) - 1
    
    for i, coeff in enumerate(coeffs):
        if coeff != 0:
            exponent = degree - i
            coeff_str = f"{int(coeff)}" if coeff.is_integer() else f"{coeff:.2f}"
            
            if coeff == 1 and exponent != 0:
                coeff_str = ""
            elif coeff == -1 and exponent != 0:
                coeff_str = "-"
            
            if exponent > 1:
                term = f"{coeff_str}x^{{{exponent}}}"
            elif exponent == 1:
                term = f"{coeff_str}x"
            else:
                term = f"{coeff_str}"
            
            terms.append(term)
    
    formatted_equation = " + ".join(terms)
    formatted_equation = formatted_equation.replace("+ -", "- ")  # Correcting signs
    
    return f"${formatted_equation}$"

def plot_from_file(filename, output_folder="images"):
    """
    Reads a file and creates a simple line plot of the data, ignoring the header.
    
    Parameters
    ----------
    filename : str
        The name of the file to read
    output_folder : str
        The folder where the plot image will be saved
    
    Returns
    -------
    None
    """
    # Create the output folder if it does not exist
    os.makedirs(output_folder, exist_ok=True)
    
    # Read the file and extract numerical data
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    # Extract header information
    method_name = lines[0].strip()
    initial_guess = float(lines[1].split(':')[1].strip())
    initial_guess_str = f"{int(initial_guess)}" if initial_guess.is_integer() else f"{initial_guess:.2f}"
    coefficients = lines[2].split('=')[1].strip()
    formatted_equation = parse_coefficients(coefficients)
    
    # Generate graph title
    grafic_title = f"{method_name} f(x) = {formatted_equation} Initial Guess = {initial_guess_str}"
    grafic_file_name = f"grafic-ex-2.2-{method_name.lower().replace(' ', '-')}-{initial_guess_str.replace('.', '')}.png"
    
    # Skip header and extract numerical values
    data = []
    for line in lines[3:]:  # Ignoring the first three lines
        parts = line.split()
        if len(parts) == 2:
            try:
                data.append([int(parts[0]), float(parts[1])])
            except ValueError:
                continue
    
    # Convert data to numpy array
    data = np.array(data)
    x = data[:, 0]
    y = data[:, 1]
    
    # Create the plot
    plt.figure(figsize=(8, 6))
    plt.plot(x, y, marker='o', linestyle='-', color='b', label='f(x) values')
    
    # Set the labels and title
    plt.xlabel("Iterations")
    plt.ylabel("f(x)")
    plt.title(grafic_title, fontsize=12)
    plt.grid()
    
    # Save the plot in the specified folder
    file_name = os.path.join(output_folder, grafic_file_name)
    plt.savefig(file_name)
    
    # Close the plot
    plt.close()

def main():
    """
    Main function to execute the plotting of data from a specified file.
    
    Returns
    -------
    None
    """
    print("Insert first data file name:")
    data_file_name = str(input().strip())

    plot_from_file(data_file_name)

    print("Insert second data file name:")
    data_file_name = str(input().strip())
    
    plot_from_file(data_file_name)

if __name__ == '__main__':
    main()
