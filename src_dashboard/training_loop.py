import torch
from synapse_utils import apply_temporal_modulation

class SynapseTrainer:
    def __init__(self, model, optimizer, loss_function, data_stream, device='cuda' if torch.cuda.is_available() else 'cpu'):
        self.model = model
        self.optimizer = optimizer
        self.loss_function = loss_function
        self.data_stream = data_stream
        self.device = device
        self.model.to(device)

    def train_cycle(self, modulation_intensity=1.0, epoch_index=0):
        self.model.train()
        cumulative_loss = 0

        for inputs, targets in self.data_stream:
            inputs, targets = inputs.to(self.device), targets.to(self.device)

            self.optimizer.zero_grad()
            predictions = self.model(inputs)
            loss = self.loss_function(predictions, targets)
            loss.backward()
            self.optimizer.step()

            apply_temporal_modulation(self.model, modulation_intensity, epoch_index)

            cumulative_loss += loss.item()

        return cumulative_loss / len(self.data_stream)

    def evaluate(self, validation_stream):
        self.model.eval()
        total_validation_loss = 0

        with torch.no_grad():
            for inputs, targets in validation_stream:
                inputs, targets = inputs.to(self.device), targets.to(self.device)
                predictions = self.model(inputs)
                loss = self.loss_function(predictions, targets)
                total_validation_loss += loss.item()

        return total_validation_loss / len(validation_stream)
