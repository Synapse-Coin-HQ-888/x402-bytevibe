import torch
import numpy as np

def generate_synaptic_waveform(length, complexity_factor=1.5):
    timeline = np.linspace(0, 1, length)
    signal = np.zeros_like(timeline)
    for i in range(1, int(length / 2)):
        signal += np.sin(2 * np.pi * (2 ** i) * timeline) / (i ** complexity_factor)
    return (signal - signal.min()) / (signal.max() - signal.min())

def apply_rhythmic_modulation(model, intensity, epoch_index):
    phase_value = generate_synaptic_waveform(1000)[epoch_index % 1000]
    for param in model.parameters():
        param.data += torch.randn_like(param.data) * intensity * phase_value

def compute_model_scale(model):
    return sum(p.numel() for p in model.parameters() if p.requires_grad)

def abbreviate_number(value):
    for suffix in ['', 'K', 'M', 'B', 'T']:
        if abs(value) < 1000.0:
            return f"{value:.1f}{suffix}"
        value /= 1000.0
    return f"{value:.1f}P"
