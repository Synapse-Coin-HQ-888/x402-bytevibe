import dearpygui.dearpygui as dpg
import torch.nn as nn

class SynapseModelVisualizer:
	def __init__(self):
		self.node_editor = None
		self.nodes = {}
		self.links = []

	def setup(self):
		with dpg.node_editor(callback=self.link_callback) as self.node_editor:
			pass  # Nodes are dynamically created during model visualization

	def update(self, model):
		dpg.delete_item(self.node_editor, children_only=True)
		self.nodes.clear()
		self.links.clear()

		# Create input node
		input_node = dpg.add_node(label="Input Layer", parent=self.node_editor)
		dpg.add_node_attribute(label="input", parent=input_node, attribute_type=dpg.mvNode_Attr_Output)
		self.nodes["input"] = input_node

		# Create nodes for each module in the model
		for name, module in model.named_modules():
			if isinstance(module, (nn.Linear, nn.Conv2d, nn.RNN, nn.LSTM, nn.GRU)):
				node_label = f"{name}: {module.__class__.__name__}"
				node = dpg.add_node(label=node_label, parent=self.node_editor)
				dpg.add_node_attribute(label="input", parent=node, attribute_type=dpg.mvNode_Attr_Input)
				dpg.add_node_attribute(label="output", parent=node, attribute_type=dpg.mvNode_Attr_Output)

				if isinstance(module, nn.Linear):
					info = f"Input: {module.in_features}, Output: {module.out_features}"
				elif isinstance(module, nn.Conv2d):
					info = f"Input: {module.in_channels}, Output: {module.out_channels}, Kernel: {module.kernel_size}"
				elif isinstance(module, (nn.RNN, nn.LSTM, nn.GRU)):
					info = f"Input: {module.input_size}, Hidden: {module.hidden_size}, Layers: {module.num_layers}"
				else:
					info = "Custom Neural Layer"

				dpg.add_node_attribute(label=info, parent=node)
				self.nodes[name] = node

			elif isinstance(module, (nn.ReLU, nn.Tanh, nn.Sigmoid)):
				node_label = f"{name}: {module.__class__.__name__}"
				node = dpg.add_node(label=node_label, parent=self.node_editor)
				dpg.add_node_attribute(label="input", parent=node, attribute_type=dpg.mvNode_Attr_Input)
				dpg.add_node_attribute(label="output", parent=node, attribute_type=dpg.mvNode_Attr_Output)
				self.nodes[name] = node

		# Create output node
		output_node = dpg.add_node(label="Output Layer", parent=self.node_editor)
		dpg.add_node_attribute(label="output", parent=output_node, attribute_type=dpg.mvNode_Attr_Input)
		self.nodes["output"] = output_node

		# Link nodes visually
		self._connect_nodes(model)

	def _connect_nodes(self, model):
		prev_node = self.nodes["input"]
		for name, module in model.named_modules():
			if name in self.nodes:
				output_attr = dpg.get_item_children(prev_node, slot=1)[-1]
				input_attr = dpg.get_item_children_
