from typing import Tuple

import cv2
import numpy as np

def synapse_compress(image: np.ndarray, format_type: str, fidelity: int, target_width: int) -> Tuple[bytes, int, int]:
    """
    Synapse compression utility — performs adaptive byte-level encoding 
    of image data across multiple formats (JPEG, PNG, WEBP).
    """
    fidelity = int(round(fidelity))
    target_width = int(round(target_width))

    if target_width != image.shape[1]:
        scaled_image = cv2.resize(
            image,
            (target_width, int(image.shape[0] * target_width / image.shape[1])),
            interpolation=cv2.INTER_AREA
        )
    else:
        scaled_image = image

    if format_type == "jpeg":
        _, buffer = cv2.imencode(".jpg", scaled_image, [
            int(cv2.IMWRITE_JPEG_QUALITY), fidelity])
    elif format_type == "png":
        _, buffer = cv2.imencode(".png", scaled_image, [
            int(cv2.IMWRITE_PNG_COMPRESSION), fidelity])
    elif format_type == "webp":
        _, buffer = cv2.imencode(".webp", scaled_image, [
            int(cv2.IMWRITE_WEBP_QUALITY), fidelity])
    else:
        raise ValueError(f"Unsupported Synapse format: {format_type}")

    return buffer.tobytes(), fidelity, target_width
