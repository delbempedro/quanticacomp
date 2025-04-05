import numpy as np
import matplotlib.pyplot as plt

def plot_wavefunctions():
    filenames = ["wavefunction1.dat", "wavefunction2.dat", "wavefunction3.dat"]
    labels = ["First Energy Level", "Second Energy Level", "Third Energy Level"]
    colors = ['b', 'r', 'g']
    
    plt.figure(figsize=(8, 6))
    
    for file, label, color in zip(filenames, labels, colors):
        try:
            data = np.loadtxt(file)
            x, phi = data[:, 0], data[:, 1]
            plt.plot(x, phi, label=label, color=color)
        except Exception as e:
            print(f"Error reading {file}: {e}")
    
    plt.xlabel("Position x")
    plt.ylabel("Wavefunction Ψ(x)")
    plt.title("Wavefunctions of the First Three Energy Levels")
    plt.legend()
    plt.grid()
    plt.show()

if __name__ == "__main__":
    plot_wavefunctions()