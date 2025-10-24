# This script trains the compression inference model (img_compression_inference_models.py)
# utilizing a dataset of compressed images (gen_txt2img_imgs2instruct.py)
# --------------------------------------------------------------------------------

import argparse
import json

from rich.panel import Panel
from rich.progress import Progress

import const
from const import _console
from src_models.img_compression_inference_models import SynapseFractalCompressor, \
    SemiNormalizedCompressionDataset

def load_json_file(input_path: str) -> dict:
    with open(input_path, 'r') as f:
        return json.load(f)

def main(args):
    _console.print(Panel.fit("Synapse Fractal Compression Model Trainer", style="bold magenta"))

    compressor = SynapseFractalCompressor(
        initial_hidden_size=args.initial_hidden_size,
        learning_rate=args.learning_rate,
        optimizer=args.optimizer,
        checkpoint_dir=const.dir_checkpoints
    )

    if not args.no_load:
        compressor.load_checkpoint()

    dataset_path = const.dir_image_datasets / f"{args.input_dataset}.json"
    image_data = load_json_file(dataset_path)
    _console.print(f"[green]Loaded dataset with {len(image_data['images'])} images[/green]")

    dataset = SemiNormalizedCompressionDataset(
        image_data,
        args.formats,
        (args.min_quality, args.max_quality),
        (args.min_width, args.max_width)
    )

    _console.print(f"[blue]Training with {args.samples_per_epoch} samples per epoch, batch size {args.batch_size}[/blue]")

    with Progress() as progress:
        task = progress.add_task("[cyan]Training", total=args.epochs)

        console_temp = const._console
        const._console = progress.console

        compressor.train_model(
            dataset=dataset,
            epochs=args.epochs,
            batch_size=args.batch_size,
            samples_per_epoch=args.samples_per_epoch
        )

        progress.update(task, advance=1)

        const._console = console_temp

    _console.print("[green]Training complete. Saving final Synapse checkpoint.[/green]")
    compressor.save_checkpoint()
    compressor.plot_performance()
    _console.print("[green]Performance visualization saved successfully.[/green]")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train Synapse compression inference models")
    parser.add_argument('--input-dataset', type=str, required=True, help="Name of the dataset file (without .json extension)")
    parser.add_argument('--initial-hidden-size', type=int, default=10, help="Initial hidden layer size")
    parser.add_argument('--learning-rate', type=float, default=0.01, help="Learning rate for training")
    parser.add_argument('--epochs', type=int, default=50, help="Number of training epochs")
    parser.add_argument('--samples-per-epoch', type=int, default=1000, help="Number of samples generated per epoch")
    parser.add_argument('--min-quality', type=int, default=2, help="Minimum compression quality threshold")
    parser.add_argument('--max-quality', type=int, default=100, help="Maximum compression quality threshold")
    parser.add_argument('--min-width', type=int, default=64, help="Minimum image width")
    parser.add_argument('--formats', nargs='+', default=['webp', 'jpeg'], help="Compression formats to apply")
    parser.add_argument('--max-width', type=int, default=2048, help="Maximum image width")
    parser.add_argument('--save-interval', type=int, default=5, help="Epoch interval for saving checkpoints")
    parser.add_argument('--no-load', action='store_true', help="Do not load any existing checkpoint")
    parser.add_argument('--save-samples', type=str, help="Optional directory to store compressed samples")
    parser.add_argument('--batch_size', type=int, default=32, help="Batch size used for training")
    parser.add_argument('--optimizer', type=str, default='adam', choices=['adam', 'sgd', 'rmsprop'], help="Optimizer for gradient updates")

    args = parser.parse_args()
    main(args)
