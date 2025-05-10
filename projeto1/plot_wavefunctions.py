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
import matplotlib.pyplot as plt

def read_wavefunctions(filename):
    """
    Reads wavefunction data from a file.

    Each wavefunction is assumed to be separated by a line with only "----------------------------------------------".
    Each data line must contain two floating-point numbers: amplitude and position.

    Parameters:
        filename (str): Path to the file containing wavefunction data.

    Returns:
        list of tuples: Each tuple contains two lists (positions, amplitudes).
    """
    with open(filename, 'r') as file:
        lines = file.readlines()

    wavefunctions = []
    positions = []
    amplitudes = []

    for line in lines:
        line = line.strip()
        if line == '----------------------------------------------':
            if positions and amplitudes:
                wavefunctions.append((positions, amplitudes))
                positions, amplitudes = [], []
        elif line:
            try:
                amplitude, position = map(float, line.split())
                positions.append(position)
                amplitudes.append(amplitude)
            except ValueError:
                continue  # Ignore improperly formatted lines

    # Append the last wavefunction if file does not end with separator
    if positions and amplitudes:
        wavefunctions.append((positions, amplitudes))

    return wavefunctions

def plot_all_wavefunctions(wavefunctions, output_file='wavefunctions.png'):
    """
    Plots all wavefunctions on the same graph and saves the plot to a PNG file.

    Parameters:
        wavefunctions (list of tuples): Each tuple contains (positions, amplitudes).
        output_file (str): Filename for the saved plot.
    """
    plt.figure(figsize=(10, 6))

    for index, (positions, amplitudes) in enumerate(wavefunctions):
        plt.plot(positions, amplitudes, label=f'Wavefunction {index + 1}')

    plt.title('All Wavefunctions')
    plt.xlabel('Position')
    plt.ylabel('Wavefunction Amplitude')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

def main():
    """
    Main function to read and plot wavefunctions from a file.
    """
    filename = 'eigenfunctions.txt'  # Adjust if your file has a different name
    wavefunctions = read_wavefunctions(filename)
    plot_all_wavefunctions(wavefunctions)

if __name__ == '__main__':
    main()

