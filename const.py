from pathlib import Path
from rich.console import Console
import matplotlib as mpl

# Define the Synapse project root directory
project_root = Path(__file__).parent

# Define and ensure all necessary directories exist
dir_instruct_datasets = project_root / Path("synapse_instruct_datasets")
dir_instruct_datasets.mkdir(exist_ok=True, parents=True)

dir_image_datasets = project_root / Path("synapse_image_datasets")
dir_image_datasets.mkdir(exist_ok=True, parents=True)

dir_audio_datasets = project_root / Path("synapse_audio_datasets")
dir_audio_datasets.mkdir(exist_ok=True, parents=True)

dir_checkpoints = project_root / Path("synapse_checkpoints")
dir_checkpoints.mkdir(exist_ok=True, parents=True)

dir_plots = project_root / Path(".synapse_plots")
dir_plots.mkdir(exist_ok=True, parents=True)

# Initialize console for formatted output
_console = Console()

# Configure Matplotlib backend to avoid PyQt dependency
mpl.use('TKAgg')  # Using TK backend to ensure compatibility and reduce external dependencies
