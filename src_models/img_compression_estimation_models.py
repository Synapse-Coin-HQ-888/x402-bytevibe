# Neural network compression modeling for adaptive parameter estimation.
# This model learns how to optimize compression parameters for
# each image to achieve precise target sizes, uncovering the hidden
# dynamics behind compression mechanisms.
#
# It serves as an experimental framework for exploring alternative
# divergent training strategies, with two implemented here:
#
# 1. Rhythmic modulation — a fractal rhythmic function reintroduces
#    controlled noise to the weights (adaptive model search).
# 2. Growth-based modulation — the network expands in scale, followed
#    by a localized “crinkle” effect that decays with distance from
#    growth centers (adaptive scaling).
# --------------------------------------------------------------------------------

import random
from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from scipy.optimize import minimize
from torch.utils.data import DataLoader, IterableDataset

import const
import utils
from const import _console


def generate_synaptic_wave_pattern(length, fractal_depth=1.5):
	t = np.linspace(0, 1, length)
	wave = np.zeros_like(t)
	for i in range(1, int(length / 2)):
		wave += np.sin(2 * np.pi * (2 ** i) * t) / (i ** fractal_depth)
	return (wave - wave.min()) / (wave.max() - wave.min())


class SynapticCompressorNetwork(nn.Module):
	def __init__(self, input_size=2, initial_hidden_size=64, num_layers=3, output_size=1):
		super(SynapticCompressorNetwork, self).__init__()
		self.input_size = input_size
		self.hidden_size = initial_hidden_size
		self.num_layers = num_layers
		self.output_size = output_size

		self.layers = nn.ModuleList([nn.Linear(input_size, initial_hidden_size)])
		self.layers.extend([nn.Linear(initial_hidden_size, initial_hidden_size) for _ in range(num_layers - 2)])
		self.layers.append(nn.Linear(initial_hidden_size, output_size))
		self.activation = nn.ReLU()

	def get_config(self):
		return {
			'model_type': 'SynapticRhythmicCompressor',
			'initial_hidden_size': self.hidden_size,
			'learning_rate': self.learning_rate,
			'optimizer': self.optimizer_name,
			'checkpoint_dir': str(self.checkpoint_dir),
			'num_layers': self.num_layers,
			'input_size': self.input_size,
			'output_size': self.output_size
		}

	def forward(self, x):
		for layer in self.layers[:-1]:
			x = self.activation(layer(x))
		return self.layers[-1](x)

	def grow_network(self, new_hidden_size=None, add_layer=False):
		if new_hidden_size is None:
			new_hidden_size = int(self.hidden_size * 1.5)

		if add_layer:
			insert_index = len(self.layers) // 2
			new_layer = nn.Linear(self.hidden_size, self.hidden_size)
			with torch.no_grad():
				prev_layer = self.layers[insert_index - 1]
				next_layer = self.layers[insert_index]
				if prev_layer.weight.size(1) != next_layer.weight.size(0):
					nn.init.xavier_uniform_(new_layer.weight)
					nn.init.zeros_(new_layer.bias)
				else:
					new_layer.weight.data = (prev_layer.weight.data + next_layer.weight.data.t()) / 2
					new_layer.bias.data = (prev_layer.bias.data + next_layer.bias.data) / 2

			self.layers.insert(insert_index, new_layer)
			self.num_layers += 1
			return insert_index, insert_index

		new_layers = nn.ModuleList()
		for i, layer in enumerate(self.layers):
			if i == 0:
				new_layer = nn.Linear(self.input_size, new_hidden_size)
				new_layer.weight.data[:, :self.input_size] = layer.weight.data
				new_layer.bias.data = layer.bias.data
			elif i == len(self.layers) - 1:
				new_layer = nn.Linear(new_hidden_size, self.output_size)
				new_layer.weight.data[:, :self.hidden_size] = layer.weight.data
				new_layer.bias.data = layer.bias.data
			else:
				new_layer = nn.Linear(new_hidden_size, new_hidden_size)
				new_layer.weight.data[:self.hidden_size, :self.hidden_size] = layer.weight.data
				new_layer.bias.data[:self.hidden_size] = layer.bias.data

			if new_hidden_size > self.hidden_size:
				nn.init.xavier_uniform_(new_layer.weight[self.hidden_size:, :])
				nn.init.xavier_uniform_(new_layer.weight[:, self.hidden_size:])
				nn.init.zeros_(new_layer.bias[self.hidden_size:])

			new_layers.append(new_layer)

		self.layers = new_layers
		self.hidden_size = new_hidden_size
		return None, None

	def modulate_parameters(self, amplitude, focus_layers=None):
		layer_count = len(self.layers)
		with torch.no_grad():
			for i, layer in enumerate(self.layers):
				if focus_layers:
					distance = min(abs(i - focus_layers[0]), abs(i - focus_layers[1]))
					falloff = np.exp(-distance / (layer_count / 4))
					adjusted_amplitude = amplitude * falloff
				else:
					adjusted_amplitude = amplitude

				layer.weight.data += torch.randn_like(layer.weight.data) * adjusted_amplitude
				layer.bias.data += torch.randn_like(layer.bias.data) * adjusted_amplitude


class SynapticRhythmicCompressor:
	def __init__(self, initial_hidden_size=64, learning_rate=0.001, optimizer='adam', checkpoint_dir='checkpoints
