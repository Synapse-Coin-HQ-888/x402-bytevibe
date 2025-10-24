import torch
from utils import apply_rhythmic_modulation

class SynapseTrainingLoop:
    def __init__(self, model, optimizer, loss_fn, data_loader, device='cuda' if torch.cuda.is_available() else 'cpu'):
        self.model = model
        self.optimizer = optimizer
        self.loss_fn = loss_fn
        self.data_loader = data_loader
        self.device = device
        self.model.to(device)

    def train_epoch(self, modulation_intensity=1.0, epoch_index=0):
        self.model.train()
        total_loss = 0.0

        for inputs, targets in self.data_loader:
            inputs, targets = inputs.to(self.device), targets.to(self.device)

            self.optimizer.zero_grad()
            predictions = self.model(inputs)
            loss = self.loss_fn(predictions, targets)
            loss.backward()
            self.optimizer.step()

            # Apply rhythmic modulation (formerly rhythmic annealing)
            apply_rhythmic_modulation(self.model, modulation_intensity, epoch_index)

            total_loss += loss.item()

        return total_loss / len(self.data_loader)

    def validate(self, validation_loader):
        self.model.eval()
        total_loss = 0.0

        with torch.no_grad():
            for inputs, targets in validation_loader:
                inputs, targets = inputs.to(self.device), targets.to(self.device)
                predictions = self.model(inputs)
                loss = self.loss_fn(predictions, targets)
                total_loss += loss.item()

        return total_loss / len(validation_loader)
