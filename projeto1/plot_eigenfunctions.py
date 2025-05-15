#!/usr/bin/env python3
"""
plot_eigenfunctions.py

Module to plot the convergence results of Matching Methods.

Dependencies:
- Use numpy to read data
- Use matplotlib to plot

Since:
- 05/2025

Authors:
- Pedro C. Delbem. <pedrodelbem@usp.br>
"""
import numpy as np
import matplotlib.pyplot as plt

def read_eigenfunctions(filename):
    states = []
    with open(filename, 'r') as f:
        lines = f.readlines()

    state_data = None
    r_vals, V_vals, psi_vals = [], [], []

    for line in lines:
        line = line.strip()
        if not line:
            if state_data is not None:
                state_data['r'] = np.array(r_vals)
                state_data['V'] = np.array(V_vals)
                state_data['psi'] = np.array(psi_vals)
                states.append(state_data)
                state_data = None
                r_vals, V_vals, psi_vals = [], [], []
            continue

        if line.startswith('State'):
            parts = line.split()
            state_num = int(parts[1])
            energy = float(parts[-1])
            state_data = {'state': state_num, 'energy': energy}
        else:
            r, V, psi = map(float, line.split())
            r_vals.append(r)
            V_vals.append(V)
            psi_vals.append(psi)

    if state_data is not None and len(r_vals) > 0:
        state_data['r'] = np.array(r_vals)
        state_data['V'] = np.array(V_vals)
        state_data['psi'] = np.array(psi_vals)
        states.append(state_data)

    return states

def normalize_potential_range(V, new_min=-25.0, new_max=25.0):
    V_min = np.min(V)
    V_max = np.max(V)
    if V_max == V_min:
        return V
    return (V - V_min) / (V_max - V_min) * (new_max - new_min) + new_min

def plot_all_together(states, psi_scale=5.0):
    plt.figure(figsize=(10, 6))

    r = states[0]['r']
    V_raw = states[0]['V']
    V_norm = normalize_potential_range(V_raw, -25, 25)

    plt.plot(r, V_norm, 'k--', label='Potencial V(r) normalizado [-25, 25]', linewidth=2)

    colors = plt.cm.viridis(np.linspace(0, 1, len(states)))

    for state, color in zip(states, colors):
        psi_scaled = state['psi'] * psi_scale + state['energy']
        plt.plot(r, psi_scaled, color=color,
                 label=f'Estado {state["state"]} (k={state["energy"]:.4f})')

    plt.xlabel('r')
    plt.ylabel('ψ(r) ampliado + energia')
    plt.title('Autofunções amplificadas e Potencial Lennard-Jones (normalizado entre -25 e 25)')
    plt.grid(True)

    # Legenda fora do gráfico (lateral direita)
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    filename = "eigenfunctions_total.txt"
    states = read_eigenfunctions(filename)
    plot_all_together(states, psi_scale=5.0)


