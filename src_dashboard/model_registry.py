from src_models.img_compression_estimation_models import FractalRhythmicCompressor

class SynapseModelHub:
    def __init__(self):
        self._registry = {
            "FractalRhythmicCompressor": FractalRhythmicCompressor
        }

    def add_model(self, identifier, model_reference):
        self._registry[identifier] = model_reference

    def load_model(self, identifier, **kwargs):
        if identifier not in self._registry:
            raise ValueError(f"Model {identifier} is not available in Synapse registry")
        return self._registry[identifier](**kwargs)

    def list_models(self):
        return list(self._registry.keys())

synapse_instance = SynapseModelHub()
