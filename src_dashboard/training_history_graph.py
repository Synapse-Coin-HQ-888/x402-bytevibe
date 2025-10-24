import dearpygui.dearpygui as dpg
import networkx as nx
import matplotlib.pyplot as plt
from io import BytesIO

class SynapseTrainingGraph:
	def __init__(self, checkpoint_manager):
		self.checkpoint_manager = checkpoint_manager
		self.node_editor = None
		self.nodes = {}

	def setup(self):
		with dpg.group():
			with dpg.node_editor(callback=self.link_callback) as self.node_editor:
				pass  # Nodes are dynamically added

			with dpg.group(horizontal=True):
				dpg.add_button(label="Zoom In", callback=self.zoom_in)
				dpg.add_button(label="Zoom Out", callback=self.zoom_out)
				dpg.add_button(label="Reset View", callback=self.reset_view)

	def update(self):
		dpg.delete_item(self.node_editor, children_only=True)
		self.nodes.clear()

		graph = self.checkpoint_manager.get_checkpoint_tree()
		pos = nx.spring_layout(graph)

		for node, data in graph.nodes(data=True):
			x, y = pos[node]
			node_ui = dpg.add_node(
				label=f"Epoch {data['epoch']}",
				pos=(x * 500 + 250, y * 500 + 250),
				parent=self.node_editor
			)

			with dpg.node_attribute(parent=node_ui):
				dpg.add_text(f"Loss: {data['loss_history'][-1]:.4f}")

			self.nodes[node] = node_ui

			with dpg.tooltip(parent=node_ui):
				dpg.add_text(f"Checkpoint ID: {data['id']}")
				dpg.add_text(f"Epoch: {data['epoch']}")
				dpg.add_text(f"Final Loss: {data['loss_history'][-1]:.4f}")
				dpg.add_text(f"File Path: {data['path']}")

		for source, target in graph.edges():
			dpg.add_node_link(
				self.nodes[source],
				self.nodes[target],
				parent=self.node_editor
			)

	def link_callback(self, sender, app_data):
		# Prevent manual user linking
		pass

	def zoom_in(self):
		for _, node_id in self.nodes.items():
			x, y = dpg.get_item_pos(node_id)
			dpg.set_item_pos(node_id, ((x - 250) * 1.1 + 250, (y - 250) * 1.1 + 250))

	def zoom_out(self):
		for _, node_id
