# This script converts an existing Synapse-generated image dataset (from gen_txt2img_imgs.py)
# into an instruction-based dataset for multimodal model training.
# It optionally applies adaptive image compression using learned estimation models
# that are included and trained within this repository.
# --------------------------------------------------------------------------------

import argparse
import base64
import json
import math
import random
import shutil
from pathlib import Path
from typing import List, Optional

import cv2
from pydantic import BaseModel, Field
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress

import const
from gen_txt2img_imgs import ImageDataset
from src_models.img_compression_estimation_models import NeuralNetworkCompressor
from utils import try_compression

# Special byte pattern marking a shift in data modality
MAGIC_BYTES = b'\xF0\x9F\x8C\x88'  # Unicode rainbow emoji marker

_console = Console()

DATASET_PREFIX = "txt2img"

USER_PROMPTS = [
    "Can you create an image of {prompt}?",
    "I'd like to see a visual of {prompt}. Can you make that?",
    "Please generate an image showing {prompt}.",
    "Could you render a depiction of {prompt}?",
    "I'm interested in {prompt}. Could you visualize that?",
    "txt2img({prompt})",
    "txt2{format}({prompt})",
]

ASSISTANT_PROMPTS = [
    "Certainly! I'll generate the requested image.",
    "Of course! Creating the image now.",
    "Sure! Generating the image for your request.",
    "Absolutely! I’ll render the image for you.",
    "I'm ready to visualize your request.",
]

out_workdir = const.dir_instruct_datasets / ".synapse_txt2img_dataset"

class CompressedImageData(BaseModel):
    file_path: str
    compression_format: str
    compression_params: dict

class InstructSample(BaseModel):
    id: int
    user_prompt: str
    image_prompt: str
    original_image_path: str
    compressed_image: Optional[CompressedImageData] = None
    seed: Optional[int] = None
    instruct_sample_path: Optional[str] = None

class InstructDataset(BaseModel):
    samples: List[InstructSample] = Field(default_factory=list)
    metadata: dict = Field(default_factory=dict)

def load_json(in_file_path: str) -> dict:
    try:
        with open(in_file_path, 'r') as f:
            return json.load(f)
    except FileNotFoundError:
        _console.print(f"[bold red]Error: File {in_file_path} not found.")
        raise
    except json.JSONDecodeError:
        _console.print(f"[bold red]Error: Invalid JSON format in {in_file_path}.")
        raise

def save_json(data: dict, out_file_path: str):
    try:
        with open(out_file_path, 'w') as f:
            json.dump(data, f, indent=2)
    except IOError:
        _console.print(f"[bold red]Error: Unable to write to {out_file_path}.")
        raise

def parse_range(value: str) -> List[int]:
    result = []
    for part in value.split(','):
        if '-' in part:
            start, end = map(int, part.split('-'))
            result.extend(range(start, end + 1))
        else:
            result.append(int(part))
    return result

def guided_compression(in_image_path,
                       out_compressed_images,
                       cfg_formats,
                       cfg_params,
                       cfg_resize,
                       byte_range,
                       console,
                       compressor):
    in_img = cv2.imread(str(in_image_path))
    quality_range = parse_range(cfg_params.get("quality", "2-100"))
    width_range = parse_range(cfg_resize) if cfg_resize else [in_img.shape[1]]
    min_bytes, max_bytes = byte_range

    max_iters = 50
    samples_per_iter = 5

    for iteration in range(max_iters):
        console.print(f"[cyan]Iteration {iteration + 1}[/cyan]")
        new_data = []

        for _ in range(samples_per_iter):
            fmt = random.choice(cfg_formats)
            quality = random.randint(min(quality_range), max(quality_range))
            width = random.randint(min(width_range), max(width_range))

            try:
                compressed, used_quality, used_width = try_compression(in_img, fmt, quality, width)
                size = len(compressed)
                new_data.append((used_quality, used_width, size))
                console.print(f"[dim]Sample: {fmt} q={used_quality} w={used_width}\t{size}b[/dim]")

                if min_bytes <= size <= max_bytes:
                    compressor.update_model(new_data)
                    return CompressedImageData(
                        file_path=save_compressed_image(in_image_path, out_compressed_images, compressed, fmt, used_quality, used_width),
                        compression_format=fmt,
                        compression_params={"quality": used_quality, "resize": used_width}
                    )

            except Exception as e:
                console.print(f"[red]Compression failed: {fmt} q={quality} w={width}\t{str(e)}[/red]")

        compressor.update_model(new_data)

        # Optimize parameters
        optimized = compressor.get_optimized_params(quality_range, width_range, (min_bytes + max_bytes) / 2)
        for quality, width in optimized:
            try:
                compressed, used_quality, used_width = try_compression(in_img, fmt, quality, width)
                size = len(compressed)
                if min_bytes <= size <= max_bytes:
                    return CompressedImageData(
                        file_path=save_compressed_image(in_image_path, out_compressed_images, compressed, fmt, used_quality, used_width),
                        compression_format=fmt,
                        compression_params={"quality": used_quality, "resize": used_width}
                    )
            except Exception as e:
                console.print(f"[red]Optimized compression failed: q={quality} w={width}\t{str(e)}[/red]")

    raise ValueError(f"Failed to compress image within {min_bytes}-{max_bytes} bytes after {max_iters} iterations")

def save_compressed_image(in_image_path: str, out_dir: Path, compressed: bytes, fmt: str, quality: int, width: int) -> str:
    stem = Path(in_image_path).stem
    out_filename = f"{stem}_{fmt}_q{quality}_w{width}"
    out_path = out_dir / f"{out_filename}.{fmt}"
    with open(out_path, "wb") as f:
        f.write(compressed)
    return str(out_path)

def create_instruct_sample(user_prompt: str, assistant_prompt: str, in_image_path: str) -> bytes:
    start = f"User: {user_prompt}\n\nAssistant: {assistant_prompt}\n\n".encode('utf-8')
    with open(in_image_path, 'rb') as f:
        image_bytes = f.read()
    end = b"\n\nHere’s the image I generated based on your request."
    return start + MAGIC_BYTES + image_bytes + MAGIC_BYTES + end

def write_instruct_sample(sample: bytes, out_path: Path):
    with open(out_path, 'wb') as f:
        f.write(sample)

def generate_dataset(src_dataset: ImageDataset, in_dataset_stem: str, cfg_formats: List[str],
                     cfg_params: dict, cfg_resize: Optional[str], byte_range: Optional[int] = None) -> InstructDataset:
    compressor = NeuralNetworkCompressor(checkpoint_dir=const.dir_checkpoints)
    compressor.load_checkpoint()

    out_dataset = InstructDataset()
    shutil.rmtree(out_workdir, ignore_errors=True)
    out_workdir.mkdir(exist_ok=True, parents=True)

    use_compression = len(cfg_formats) > 0 and cfg_formats != ["none"]
    out_compressed_images = out_workdir / "compressed_images" if use_compression else None
    if use_compression:
        out_compressed_images.mkdir(exist_ok=True)

    out_instruct_samples = const.dir_instruct_datasets / f"{DATASET_PREFIX}_{in_dataset_stem}"
    out_instruct_samples.mkdir(exist_ok=True)

    with Progress() as progress:
        task = progress.add_task("[cyan]Generating instruction samples...", total=len(src_dataset.images))
        for i, entry in enumerate(src_dataset.images):
            img_prompt = entry.prompt
            user_prompt = random.choice(USER_PROMPTS).replace("{prompt}", img_prompt.lower())
            assistant_prompt = random.choice(ASSISTANT_PROMPTS)

            in_img_path = const.dir_image_datasets / Path(entry.image_path)
            if not in_img_path.exists():
                progress.console.print(f"[yellow]Skipping missing image: {entry.image_path}")
                continue

            compressed_image = None
            if use_compression:
                compressed_image = guided_compression(
                    in_image_path=in_img_path,
                    out_compressed_images=out_compressed_images,
                    cfg_formats=cfg_formats,
                    cfg_params=cfg_params,
                    cfg_resize=cfg_resize,
                    byte_range=byte_range,
                    console=progress.console,
                    compressor=compressor
                )

            used_image_path = compressed_image.file_path if compressed_image else in_img_path
            user_prompt = user_prompt.replace("{compression}", compressed_image.compression_format if use_compression else "img")

            instruct_sample = create_instruct_sample(user_prompt, assistant_prompt, used_image_path)
            out_sample_path = out_instruct_samples / f"sample_{i:04d}.bin"
            write_instruct_sample(instruct_sample, out_sample_path)

            sample_entry = InstructSample(
                id=i,
                user_prompt=user_prompt,
                image_prompt=img_prompt,
                original_image_path=entry.image_path,
                compressed_image=compressed_image,
                seed=entry.seed,
                instruct_sample_path=str(out_sample_path.relative_to(const.dir_instruct_datasets))
            )
            out_dataset.samples.append(sample_entry)

            progress.update(task, advance=1, description=f"[cyan]Processing {out_sample_path.name}")
            progress.console.print(f"[green]{out_sample_path.name}[/green] -> [yellow]{entry.prompt}[/yellow]")

    compressor.save_checkpoint()
    return out_dataset

def main(args):
    _console.print(Panel.fit("Synapse Instructional Image Dataset Generator", style="bold magenta"))
    use_compression = args.compression != "none"

    try:
        in_stem = Path(args.input_dataset).stem
        in_json = const.dir_image_datasets / f"{in_stem}.json"
        out_json = const.dir_instruct_datasets / f"txt2img_{in_stem}.json"

        _console.log(f"[green]Loading dataset from {in_json}...")
        src_dataset = ImageDataset(**load_json(in_json.as_posix()))
        _console.log(f"[green]Loaded {len(src_dataset.images)} images")

        cfg_params = {"quality": args.quality, "compression_level": args.compression_level}
        out_dataset = generate_dataset(src_dataset, in_stem, args.compression, cfg_params, args.resize, (args.min_bytes, args.max_bytes))

        out_dataset.metadata = {
            "num_samples": len(out_dataset.samples),
            "magic_bytes": base64.b64encode(MAGIC_BYTES).decode('ascii'),
            "format_version": "1.0",
            "description": "Instructional dataset linking text prompts with compressed images (Synapse version)",
            "compression": args.compression if use_compression else "none",
            "compression_params": cfg_params if use_compression else {},
            "use_compression": use_compression,
            "resize": args.resize
        }

        save_json(out_dataset.dict(), out_json)
        _console.print(Panel.fit("Dataset successfully generated!", style="bold green"))
        _console.print(f"[yellow]Samples: {len(out_dataset.samples)}")
        _console.print(f"[yellow]Output file: {out_json}")
        _console.print(f"[yellow]Working directory: {out_workdir}")

    except Exception as e:
        import traceback
        _console.print(f"[bold red]Error: {e}")
        _console.print(Panel(traceback.format_exc(), title="Error Details", border_style="red"))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate an instructional dataset with optional image compression (Synapse)")
    parser.add_argument("--input-dataset", type=str, required=True, help="Input dataset filename (without .json)")
    parser.add_argument("--output-dir", type=str, default=const.dir_instruct_datasets.relative_to(const.project_root), help="Output directory path")
    parser.add_argument("--compression", type=str, nargs='+', choices=["jpeg", "png", "webp"], default=["none"], help="Compression format(s)")
    parser.add_argument("--quality", type=str, default="90", help="JPEG/WebP quality (0–100, supports ranges/commas)")
    parser.add_argument("--compression-level", type=str, default="3", help="PNG compression level (0–9, supports ranges/commas)")
    parser.add_argument("--resize", type=str, help="Resize widths (supports ranges/commas)")
    parser.add_argument("--min-bytes", type=int, default=0, help="Minimum compressed image size (bytes)")
    parser.add_argument("--max-bytes", type=int, default=math.inf, help="Maximum compressed image size (bytes)")

    args = parser.parse_args()

    if args.output_dir:
        const.dir_instruct_datasets = Path(args.output_dir)
        const.dir_instruct_datasets.mkdir(exist_ok=True, parents=True)

    main(args)
