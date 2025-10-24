import torch
from torch.utils.data import Dataset, DataLoader
import json
import cv2
import random
from pathlib import Path
import const
import utils

class SynapseCompressionDataset(Dataset):
    def __init__(self, data_file, encodings, fidelity_range, resolution_range):
        self.data = self._load_json(data_file)
        self.encodings = encodings
        self.fidelity_range = fidelity_range
        self.resolution_range = resolution_range
        self.max_fidelity = max(fidelity_range)

    def _load_json(self, file_path):
        with open(file_path, "r") as f:
            return json.load(f)

    def __len__(self):
        return len(self.data["images"])

    def __getitem__(self, index):
        entry = self.data["images"][index]
        img_path = const.dir_image_datasets / entry["image_path"]
        img = cv2.imread(str(img_path))

        if img is None:
            return self.__getitem__(random.randint(0, len(self) - 1))

        encoding = random.choice(self.encodings)
        fidelity = random.randint(*self.fidelity_range)
        resolution = random.randint(*self.resolution_range)

        try:
            compressed_data = utils.try_compression(img, encoding, fidelity, resolution)
            size_bytes = len(compressed_data)

            normalized_fidelity = fidelity / self.max_fidelity

            return (
                torch.tensor([normalized_fidelity, resolution], dtype=torch.float32),
                torch.tensor([size_bytes], dtype=torch.float32),
            )
        except Exception:
            return self.__getitem__(random.randint(0, len(self) - 1))


def synapse_data_loader(data_file, encodings, fidelity_range, resolution_range, batch_size=32):
    dataset = SynapseCompressionDataset(data_file, encodings, fidelity_range, resolution_range)
    return DataLoader(dataset, batch_size=batch_size, shuffle=True)
