import numpy as np
import matplotlib.pyplot as plt

def generate_synapse_waveform(size, complexity_factor=1.5):
    """
    Generates a fractal-inspired rhythmic waveform — a Synapse pattern generator
    emulating recursive harmonic structures with tunable complexity.
    """
    timeline = np.linspace(0, 1, size)
    waveform = np.zeros_like(timeline)
    for i in range(1, int(size / 2)):
        waveform += np.sin(2 * np.pi * (2 ** i) * timeline) / (i ** complexity_factor)
    return (waveform - waveform.min()) / (waveform.max() - waveform.min())

# Parameters
size = 1000
complexity_factor = 1.5

# Generate the Synapse rhythmic waveform
synapse_wave = generate_synapse_waveform(50, complexity_factor)

# Plot configuration
plt.figure(figsize=(12, 6))
plt.plot(np.linspace(0, 1, size), synapse_wave, color='black')
plt.title(f'Synapse Waveform (Complexity Factor: {complexity_factor})')
plt.xlabel('Temporal Axis')
plt.ylabel('Amplitude')
plt.grid(True, linestyle='--', alpha=0.7)

# Add a mid-reference line at y = 0.5
plt.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)

# Render the visualization
plt.tight_layout()
plt.show()
