"""
# SynapseNet: Byte-Level Hybrid of Phi-3.5 and RetNet

SynapseNet is an experimental framework that reimagines Phi-3.5 as a byte-level model, infused with RetNet-style architectural principles. This design allows it to process raw byte inputs, enabling compatibility with virtually any file type.

## Core Capabilities:

1. **Byte-Level Modeling**: Directly handles raw byte sequences for universal file-type comprehension.
2. **RetNet Mechanisms**: Integrates RetNet-inspired multi-scale decay and normalization methods for efficient long-term context modeling.
3. **Dual Execution Modes**: Features both parallel training and recurrent inference modes.
4. **Targeted Fine-Tuning**: Allows selective updates (e.g., embedding and post-attention normalization) while maintaining most pretrained Phi-3.5 weights.
5. **DoRA Adaptation (Decomposed Low-Rank)**: Adds efficient parameter adaptation to self-attention output layers, optimizing new learning without overwriting existing representations.

## Implementation Overview:

- **Pretrained Weight Utilization**: Reuses frozen parameters from the original Phi-3.5 network.
- **Configurable DoRA Layers**: Allows precise selection of modules to receive DoRA-based adaptation.
- **Flexible Structure**: Supports both retention and attention mechanisms as interchangeable modes.
- **Optional Untied Embeddings**: Enables distinct input/output embeddings for flexible encoding schemes.

## Training and Inference:

- Provides customizable training loops with dynamic learning rate scheduling.
- Supports full training or checkpoint-based fine-tuning.
- Includes text generation for completion and evaluation.

## Research Goals:

- Examine the impact of retention-based modeling in byte-level architectures.
- Optimize efficiency by using dual-mode computation.
- Explore general-purpose byte-level modeling for diverse data formats.

> **Note:** SynapseNet is an experimental research prototype intended for exploration. It demonstrates combining pretrained large models with emerging hybrid architectures and lightweight adaptation strategies.

Author: Josef Albers  
Date: Aug 28, 2024
"""

import glob
import json
import math
import time
from datetime import datetime
from types import SimpleNamespace

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from huggingface_hub import snapshot_download
from datasets import load_dataset
from lion_pytorch import Lion  # Ensure you have installed lion-pytorch
from torch import autocast

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

class Tokenizer:
	def __init__(self, file_path=None):
		if file_path is None:
			self.vocab = list(range(256))
		else:
			with open(file_path, 'r') as f:
				content = f.read().lower().encode('utf-8')
			self.vocab = sorted(set(content))
		self.vocab_size = len(self.vocab)
		self.byte_to_index = {byte: index for index, byte in enumerate(self.vocab)}
		self.index_to_byte = {index: byte for index, byte in enumerate(self.vocab)}

	def encode(self, text):
		byte_seq = text.encode('utf-8')
		return [self.byte_to_index.get(byte, 0) for byte in byte_seq]

	def decode(self, indices):
		byte_seq = bytes([self.index_to_byte.get(index, 0) for index in indices])
		return byte_seq.decode('utf-8', errors='ignore')

# Other model components (SuRoPE, Phi3Attention, Phi3Retention, Phi3MLP, Phi3DecoderLayer, Phi3Model, Phi3ForCausalLM, DoRALinear)
# remain identical in structure — only their docstrings or top-level descriptions have been updated for the Synapse rebrand.

# Usage Example:
# python synapsenet.py

if __name__ == "__main__":
	main()
