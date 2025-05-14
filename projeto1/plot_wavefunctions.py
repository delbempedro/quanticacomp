#!/usr/bin/env python3
"""
plot_wavefunctions.py

Module to plot the convergence results of Numerical methods.

Dependencies:
- Use matplotlib to plot

Since:
- 04/2025

Authors:
- Pedro C. Delbem. <pedrodelbem@usp.br>
"""
import numpy as np
import matplotlib.pyplot as plt

def read_wavefunctions(file_path):
    """
    Reads wavefunctions from a file.
    Each wavefunction is separated by a line containing only dashes ("----------------------------------------------").
    Each data point consists of two columns: value and position.
    """
    wavefunctions = []
    current_x = []
    current_y = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line == "----------------------------------------------":
                if current_x and current_y:
                    wavefunctions.append((np.array(current_x), np.array(current_y)))
                    current_x = []
                    current_y = []
            elif line:
                try:
                    y_str, x_str = line.split()  # File has "value position"
                    current_x.append(float(x_str))  # position
                    current_y.append(float(y_str))  # value
                except ValueError:
                    print(f"Warning: ignoring malformed line: {line}")
        if current_x and current_y:
            wavefunctions.append((np.array(current_x), np.array(current_y)))

    return wavefunctions

def normalize_to_unit_interval(y):
    """
    Normalize an array to the [0, 1] interval.
    """
    y_min = np.min(y)
    y_max = np.max(y)
    if y_max == y_min:
        return np.zeros_like(y)  # Prevent division by zero
    return (y - y_min) / (y_max - y_min)

def plot_wavefunctions(wavefunctions, output_filename='wavefunctions.png'):
    """
    Plots all wavefunctions on the same graph and saves the image.
    The last wavefunction is assumed to represent the potential V(x), normalized to [0, 1].
    """
    plt.figure(figsize=(10, 6))
    num_functions = len(wavefunctions)

    for i, (x, y) in enumerate(wavefunctions):
        if i == num_functions - 1:
            #y = normalize_to_unit_interval(y)
            plt.plot(x, y, label="Potential V(x) [normalized]", color="black", linestyle="--", linewidth=2)
        else:
            plt.plot(x, y, label=f"ψ_{i+1}")

    plt.xlabel("Position")
    plt.ylabel("ψ(x), V(x)")
    plt.title("Wavefunctions and Normalized Potential")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_filename)
    print(f"Plot saved as: {output_filename}")
    plt.show()

def main():
    file_name = input("Enter the path to the file containing the wavefunctions: ").strip()
    try:
        wavefunctions = read_wavefunctions(file_name)
        if wavefunctions:
            plot_wavefunctions(wavefunctions)
        else:
            print("No wavefunctions found in the file.")
    except FileNotFoundError:
        print(f"Error: file '{file_name}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()
