from pathlib import Path
from typing import List

from src_synapse.checkpoint_controller import CheckpointController
from src_synapse.training_progress_plot import TrainingProgressPlot

class SynapseRun:
    def __init__(self, session_name, model_label, directory):
        self.session_name = session_name
        self.session_path = None
        self.model_label = model_label

        self.model_instance = None
        self.optimizer_instance = None
        self.loss_function = None
        self.data_stream = None
        self.training_cycle = None
        self.epoch_index = 0
        self.total_epochs = 0
        self.starting_epoch = 0
        self.is_training = False
        self.model_display = None
        self.training_progress_plot = None
        self.annealing_intensity = 1.0
        self.loss_log = []

        if directory:
            self.initialize()

    def initialize(self):
        self.session_path = Path(self.session_name)
        self.checkpoint_controller = CheckpointController(self.session_path / "checkpoints")


class SynapseRunManager:
    def __init__(self, sessions_dir="sessions"):
        self.sessions_dir = Path(sessions_dir)
        self.sessions_dir.mkdir(parents=True, exist_ok=True)
        self.sessions: List[SynapseRun] = self.load_sessions()

    def load_sessions(self) -> List[SynapseRun]:
        sessions = []
        for session_path in self.sessions_dir.iterdir():
            if session_path.is_dir():
                SynapseRun(session_path.name, None, session_path)
        return sessions

    def open_session(self, session_name):
        if session_name not in self.sessions:
            raise ValueError(f"Session '{session_name}' is not found")

    def create_session(self, session_name, model_label):
        if session_name in self.sessions:
            raise ValueError(f"Session '{session_name}' already exists")

        session_path = self.sessions_dir / session_name
        session_path.mkdir(parents=True, exist_ok=True)

        session = SynapseRun(session_name, model_label, None)
        self.sessions[session_name] = session

    def remove_session(self, session_name):
        if session_name not in self.sessions:
            raise ValueError(f"Session '{session_name}' does not exist")

        session_path = self.sessions_dir / session_name
        for entry in session_path.iterdir():
            if entry.is_file():
                entry.unlink()
            elif entry.is_dir():
                for subentry in entry.iterdir():
                    subentry.unlink()
                entry.rmdir()
        session_path.rmdir()

        del self.sessions[session_name]

synapse_manager = SynapseRunManager()
