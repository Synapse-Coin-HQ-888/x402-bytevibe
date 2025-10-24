import os
import json
import torch
from pathlib import Path
import networkx as nx

class SynapseCheckpointManager:
	def __init__(self, directory="checkpoints"):
		self.directory = Path(directory)
		self.directory.mkdir(parents=True, exist_ok=True)
		self.history_path = self.directory / "history.json"
		self.graph = nx.DiGraph()
		self._load_history()

	def _load_history(self):
		if self.history_path.exists():
			with open(self.history_path, "r") as file:
				history_data = json.load(file)
				for record in history_data:
					self.graph.add_node(record['id'], **record)
					if record['parent_id'] is not None:
						self.graph.add_edge(record['parent_id'], record['id'])
		else:
			self.graph.clear()

	def _save_history(self):
		history_data = [self.graph.nodes[node] for node in self.graph.nodes()]
		with open(self.history_path, "w") as file:
			json.dump(history_data, file, indent=2)

	def save(self, model, optimizer, epoch, loss_log, parent_id=None):
		checkpoint_id = self.graph.number_of_nodes()
		checkpoint_file = self.directory / f"synapse_ckpt_{checkpoint_id}.pth"

		torch.save({
			"model_state_dict": model.state_dict(),
			"optimizer_state_dict": optimizer.state_dict(),
			"epoch": epoch,
			"loss_log": loss_log
		}, checkpoint_file)

		checkpoint_meta = {
			"id": checkpoint_id,
			"epoch": epoch,
			"loss_log": loss_log,
			"parent_id": parent_id,
			"path": str(checkpoint_file)
		}

		self.graph.add_node(checkpoint_id, **checkpoint_meta)
		if parent_id is not None:
			self.graph.add_edge(parent_id, checkpoint_id)

		self._save_history()
		return checkpoint_id

	def load(self, checkpoint_id):
		if checkpoint_id not in self.graph.nodes:
			raise ValueError(f"Checkpoint with id {checkpoint_id} not found")

		checkpoint_meta = self.graph.nodes[checkpoint_id]
		checkpoint_data = torch.load(checkpoint_meta["path"])
		return checkpoint_data, checkpoint_meta

	def list_all(self):
		return [self.graph.nodes[node] for node in self.graph.nodes()]

	def get_tree(self):
		return self.graph

	# Optional extensions
	# def compare(self, id1, id2):
	#     c1, c2 = self.graph.nodes[id1], self.graph.nodes[id2]
	#     return {
	#         "epoch_diff": c2["epoch"] - c1["epoch"],
	#         "loss_diff": c2["loss_log"][-1] - c1["loss_log"][-1],
	#         "path_1": c1["path"],
	#         "path_2": c2["path"]
	#     }

	# def best_checkpoint(self, metric="loss"):
	#     if metric == "loss":
	#         best_id = min(self.graph.nodes, key=lambda x: self.graph.nodes[x]["loss_log"][-1])
	#     else:
	#         raise ValueError(f"Unsupported metric: {metric}")
	#     return self.graph.nodes[best_id]

	# def prune(self, keep_top_n=5, metric="loss"):
	#     sorted_nodes = sorted(self.graph.nodes, key=lambda x: self.graph.nodes[x]["loss_log"][-1])
	#     to_remove = sorted_nodes[keep_top_n:]
	#     for cid in to_remove:
	#         ckpt_path = Path(self.graph.nodes[cid]["path"])
	#         if ckpt_path.exists():
	#             ckpt_path.unlink()
	#         self.graph.remove_node(cid)
	#     self._save_history()
