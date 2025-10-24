"""
# RetNPhi: Bytewise Fusion of Phi-3.5 and RetNet

RetNPhi is a research prototype that converts Phi-3.5 into a byte-native language model while weaving in RetNet-like ideas. By operating directly on raw bytes, the system can ingest virtually any file format without a tokenizer tailored to text.

## Highlights

1. **Byte-Level I/O**: Consumes raw byte streams, enabling broad applicability across arbitrary file types.
2. **RetNet-Inspired Design**: Integrates multi-scale exponential decay and group normalization to efficiently capture long-range dependencies.
3. **Dual Execution Modes**: Parallel mode for fast training and recurrent mode for streamlined inference.
4. **Targeted Fine-Tuning**: Updates only selected components (e.g., token embedding and post-attention norms), leaving most Phi-3.5 weights frozen.
5. **DoRA (Weight-Decomposed LoRA)**: Applies DoRA to self-attention output projections for compact adaptation without discarding pretrained knowledge.

## Implementation Notes

- **Frozen Backbone**: Reuses the original Phi-3.5 weights for the majority of layers.
- **Configurable DoRA**: Precisely choose which layers/targets receive DoRA adapters.
- **Switchable Mechanisms**: Toggle between retention-based and standard attention paths.
- **Optional Untied Embeddings**: Use distinct input/output embeddings when desired.

## Training & Inference

- Efficient training loop with configurable learning rate schedules.
- Supports training from scratch or continued fine-tuning from checkpoints.
- Includes a simple generation utility for completion tasks.

## Objectives

- Probe retention-style mechanisms in a byte-focused Phi variant.
- Exploit dual-mode processing to balance training throughput and inference efficiency.
- Move toward a universal byte-model able to handle arbitrary file types.

> This is an exploratory codebase intended for experimentation, not production deployment. It showcases how pretrained models can be combined with novel architectural choices and parameter-efficient tuning techniques.

Author: Josef Albers
Date: Aug 28, 2024
"""

import glob
import json
import math
import time
from datetime import datetime
from types import SimpleNamespace

import fire
import mlx.core as mx
import mlx.nn as nn
import mlx.optimizers as optim
import numpy as np
from huggingface_hub import snapshot_download
from mlx.utils import tree_flatten, tree_unflatten

from datasets import load_dataset


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
        return [self.byte_to_index[byte] for byte in byte_seq]

    def decode(self, indices):
        byte_seq = bytes(self.index_to_byte[index] for index in indices)
        return byte_seq.decode('utf-8', errors='ignore')


class SuRoPE(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.dim = config.hidden_size // config.num_attention_heads
        self.original_max_position_embeddings = config.original_max_position_embeddings
        self.rope_theta = config.rope_theta
        self.scaling_factor = math.sqrt(
            1
            + math.log(
                config.max_position_embeddings / config.original_max_position_embeddings
            )
            / math.log(config.original_max_position_embeddings)
        )
        self._long_factor = mx.array(config.rope_scaling["long_factor"], dtype=mx.float32)
        self._short_factor = mx.array(config.rope_scaling["short_factor"], dtype=mx.float32)

    def __call__(self, q, k, position_ids):
        cos, sin = self._get_cos_sin(position_ids)
        q = (q * cos) + (self._rotate_half(q) * sin)
        k = (k * cos) + (self._rotate_half(k) * sin)
        return q, k

    def _get_cos_sin(self, position_ids):
        su_factor = self._short_factor
        position_ids_expanded = position_ids[:, None, :]
        inv_freq = 1.0 / (
            su_factor * self.rope_theta ** (mx.arange(0, self.dim, 2, dtype=mx.float32) / self.dim)
        )
        inv_freq_expanded = mx.repeat(inv_freq[None, :, None], position_ids.shape[0], axis=0)
        freqs = (inv_freq_expanded @ position_ids_expanded).transpose(0, 2, 1)
        emb = mx.concatenate([freqs, freqs], axis=-1)
        cos = mx.expand_dims(mx.cos(emb) * self.scaling_factor, axis=1)
        sin = mx.expand_dims(mx.sin(emb) * self.scaling_factor, axis=1)
        return cos, sin

    def _rotate_half(self, x):
        midpoint = x.shape[-1] // 2
        x1, x2 = x[..., :midpoint], x[..., midpoint:]
        return mx.concatenate([-x2, x1], axis=-1)


class Phi3Attention(nn.Module):
    def __init__(self, config):
        super().__init__()
        dim = config.hidden_size
        self.n_heads = n_heads = config.num_attention_heads
        self.n_kv_heads = n_kv_heads = config.num_key_value_heads
        self.num_hidden_layers = config.num_hidden_layers
        self.head_dim = head_dim = config.hidden_size // n_heads
        self.scale = head_dim ** -0.5
        chop_1 = self.n_heads * self.head_dim
        chop_2 = chop_1 + self.n_kv_heads * self.head_dim
        self.chop = [chop_1, chop_2]
        op_size = n_heads * head_dim + 2 * (n_kv_heads * head_dim)
        self.qkv_proj = nn.Linear(dim, op_size, bias=False)
        self.o_proj = nn.Linear(n_heads * head_dim, dim, bias=False)
        self.rope = SuRoPE(config)

    def __call__(self, x, position_ids, attention_mask, cache, use_recurrent_mode):
        B, L, _ = x.shape
        qkv = self.qkv_proj(x)
        q, k, v = mx.split(qkv, self.chop, axis=-1)
        q = q.reshape(B, L, self.n_heads, -1).transpose(0, 2, 1, 3)
        k = k.reshape(B, L, self.n_kv_heads, -1).transpose(0, 2, 1, 3)
        v = v.reshape(B, L, self.n_kv_heads, -1).transpose(0, 2, 1, 3)
        if cache is None:
            position_ids = mx.arange(q.shape[2], dtype=mx.float32)[None] if position_ids is None else position_ids
            q, k = self.rope(q, k, position_ids)
            mask = mx.triu(mx.full((v.shape[2], v.shape[2]), -mx.inf), k=1)
            if attention_mask is not None:
                mask += mx.where(attention_mask[:, :, None] * attention_mask[:, None, :] == 1, 0, -mx.inf)
                mask = mx.expand_dims(mask, 1)
            else:
                mask = mask[None, None]
        else:
            past_k, past_v, past_p, past_m = cache
            position_ids = past_p[:, -1:] + 1
            mask = mx.pad(past_m[:, :, -1:, :], ((0, 0), (0, 0), (0, 0), (0, 1)))
            q, k = self.rope(q, k, position_ids)
            k = mx.concatenate([past_k, k], axis=2)
            v = mx.concatenate([past_v, v], axis=2)
        cache = (k, v, position_ids, mask)
        w = (q * self.scale) @ k.transpose(0, 1, 3, 2)
        w += mask
        w = mx.softmax(w, axis=-1)
        o = w @ v
        o = o.transpose(0, 2, 1, 3).reshape(B, L, -1)
        return self.o_proj(o).astype(x.dtype), cache


class Phi3Retention(nn.Module):
    def __init__(self, config):
        super().__init__()
        self.dim = dim = config.hidden_size
        self.n_heads = n_heads = config.num_attention_heads
        self.n_kv_heads = n_kv_heads = config.num_key_value_heads
        self.num_hidden_layers = config.num_hidden_layers
        self.head_dim = head_dim = config.hidden_size // n_heads
        self.scale = head_dim ** -0.5
        chop_1 = self.n_heads * self.head_dim
        chop_2 = chop_1 + self.n_kv_heads * self.head_dim
        self.chop = [chop_1, chop_2]
        op_size = n_he
