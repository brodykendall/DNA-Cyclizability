from imports import *
from torch_geometric.explain import GNNExplainer
from torch_geometric.explain import Explanation
from torch_geometric.explain.algorithm.utils import clear_masks, set_masks
from torch_geometric.explain.config import MaskType, ModelMode, ModelTaskLevel
from torch.nn.parameter import Parameter

from math import sqrt


class CustomExplainer(GNNExplainer):
    def __init__(self, lr):
        super(CustomExplainer, self).__init__()
        self.lr = lr

    def forward(
        self,
        model: torch.nn.Module,
        data,
        *,
        target: torch.Tensor,
        index = None,
        **kwargs,
    ):
        self._train(model, data, target=target, index=index, **kwargs)

        node_mask = self._post_process_mask(
            self.node_mask,
            self.hard_node_mask,
            apply_sigmoid=True,
        )
        edge_mask = self._post_process_mask(
            self.edge_mask,
            self.hard_edge_mask,
            apply_sigmoid=True,
        )

        self._clean_model(model)

        return Explanation(node_mask=node_mask, edge_mask=edge_mask, edge_index=data.edge_index)
    
    def _train(
        self,
        model: torch.nn.Module,
        data,
        *,
        target: torch.Tensor,
        index = None,
        **kwargs,
    ):
        x, edge_index = data.x, data.edge_index
        self._initialize_masks(x, edge_index)

        parameters = []
        if self.node_mask is not None:
            parameters.append(self.node_mask)
        if self.edge_mask is not None:
            set_masks(model, self.edge_mask, edge_index, apply_sigmoid=True)
            parameters.append(self.edge_mask)

        optimizer = torch.optim.Adam(parameters, lr=self.lr)

        for i in range(self.epochs):
            optimizer.zero_grad()

            h = x if self.node_mask is None else x * self.node_mask.sigmoid()
            cur_data = data
            cur_data.x = h
            y_hat, y = model(cur_data, **kwargs), target

            if index is not None:
                y_hat, y = y_hat[index], y[index]

            loss = self._loss(y_hat, y)

            loss.backward()
            optimizer.step()

            # In the first iteration, we collect the nodes and edges that are
            # involved into making the prediction. These are all the nodes and
            # edges with gradient != 0 (without regularization applied).
            if i == 0 and self.node_mask is not None:
                if self.node_mask.grad is None:
                    raise ValueError("Could not compute gradients for node "
                                     "features. Please make sure that node "
                                     "features are used inside the model or "
                                     "disable it via `node_mask_type=None`.")
                self.hard_node_mask = self.node_mask.grad != 0.0
            if i == 0 and self.edge_mask is not None:
                if self.edge_mask.grad is None:
                    raise ValueError("Could not compute gradients for edges. "
                                     "Please make sure that edges are used "
                                     "via message passing inside the model or "
                                     "disable it via `edge_mask_type=None`.")
                self.hard_edge_mask = self.edge_mask.grad != 0.0

    # def _initialize_masks(self, data, model):
    #     node_mask_type = self.explainer_config.node_mask_type
    #     edge_mask_type = self.explainer_config.edge_mask_type

    #     device = data.x.device
    #     # (N, F), E = x.size(), edge_index.size(1)
    #     (N, F), E = data.x.size(), data.edge_index.size(1)
    #     # E = 0
    #     # for e_l in model.edge_lengths:
    #     #     e_l_attr_name = f"edge_index_{e_l}"
    #     #     E += getattr(data, e_l_attr_name, None).size(1)

    #     std = 0.1
    #     if node_mask_type is None:
    #         self.node_mask = None
    #     elif node_mask_type == MaskType.object:
    #         self.node_mask = Parameter(torch.randn(N, 1, device=device) * std)
    #     elif node_mask_type == MaskType.attributes:
    #         self.node_mask = Parameter(torch.randn(N, F, device=device) * std)
    #     elif node_mask_type == MaskType.common_attributes:
    #         self.node_mask = Parameter(torch.randn(1, F, device=device) * std)
    #     else:
    #         assert False

    #     if edge_mask_type is None:
    #         self.edge_mask = None
    #     elif edge_mask_type == MaskType.object:
    #         std = torch.nn.init.calculate_gain('relu') * sqrt(2.0 / (2 * N))
    #         self.edge_mask = Parameter(torch.randn(E, device=device) * std)
    #     else:
    #         assert False


class CustomExplainer2(GNNExplainer):
    def __init__(self, lr):
        super(CustomExplainer2, self).__init__()
        self.lr = lr

    def forward(
        self,
        model: torch.nn.Module,
        data,
        *,
        target: torch.Tensor,
        index = None,
        **kwargs,
    ):
        self._train(model, data, target=target, index=index, **kwargs)

        node_mask = self._post_process_mask(
            self.node_mask,
            self.hard_node_mask,
            apply_sigmoid=True,
        )
        edge_mask = self._post_process_mask(
            self.edge_mask,
            self.hard_edge_mask,
            apply_sigmoid=True,
        )

        self._clean_model(model)

        return Explanation(node_mask=node_mask, edge_mask=edge_mask, edge_index=data.edge_index)
    
    def _train(
        self,
        model: torch.nn.Module,
        data,
        *,
        target: torch.Tensor,
        index = None,
        **kwargs,
    ):
        x, edge_index = data.x, data.edge_index
        self._initialize_masks(x, edge_index)

        parameters = []
        if self.node_mask is not None:
            parameters.append(self.node_mask)
        if self.edge_mask is not None:
            set_masks(model, self.edge_mask, edge_index, apply_sigmoid=True)
            parameters.append(self.edge_mask)

        optimizer = torch.optim.Adam(parameters, lr=self.lr)

        for i in range(self.epochs):
            optimizer.zero_grad()

            h = x if self.node_mask is None else x * self.node_mask.sigmoid()
            # cur_data = data
            # cur_data.x = h
            y_hat, y = model(h, **kwargs), target

            if index is not None:
                y_hat, y = y_hat[index], y[index]

            loss = self._loss(y_hat, y)

            loss.backward(retain_graph=True)
            optimizer.step()

            # In the first iteration, we collect the nodes and edges that are
            # involved into making the prediction. These are all the nodes and
            # edges with gradient != 0 (without regularization applied).
            if i == 0 and self.node_mask is not None:
                if self.node_mask.grad is None:
                    raise ValueError("Could not compute gradients for node "
                                     "features. Please make sure that node "
                                     "features are used inside the model or "
                                     "disable it via `node_mask_type=None`.")
                self.hard_node_mask = self.node_mask.grad != 0.0
            if i == 0 and self.edge_mask is not None:
                if self.edge_mask.grad is None:
                    raise ValueError("Could not compute gradients for edges. "
                                     "Please make sure that edges are used "
                                     "via message passing inside the model or "
                                     "disable it via `edge_mask_type=None`.")
                self.hard_edge_mask = self.edge_mask.grad != 0.0

    # def _initialize_masks(self, data, model):
    #     node_mask_type = self.explainer_config.node_mask_type
    #     edge_mask_type = self.explainer_config.edge_mask_type

    #     device = data.x.device
    #     # (N, F), E = x.size(), edge_index.size(1)
    #     (N, F), E = data.x.size(), data.edge_index.size(1)
    #     # E = 0
    #     # for e_l in model.edge_lengths:
    #     #     e_l_attr_name = f"edge_index_{e_l}"
    #     #     E += getattr(data, e_l_attr_name, None).size(1)

    #     std = 0.1
    #     if node_mask_type is None:
    #         self.node_mask = None
    #     elif node_mask_type == MaskType.object:
    #         self.node_mask = Parameter(torch.randn(N, 1, device=device) * std)
    #     elif node_mask_type == MaskType.attributes:
    #         self.node_mask = Parameter(torch.randn(N, F, device=device) * std)
    #     elif node_mask_type == MaskType.common_attributes:
    #         self.node_mask = Parameter(torch.randn(1, F, device=device) * std)
    #     else:
    #         assert False

    #     if edge_mask_type is None:
    #         self.edge_mask = None
    #     elif edge_mask_type == MaskType.object:
    #         std = torch.nn.init.calculate_gain('relu') * sqrt(2.0 / (2 * N))
    #         self.edge_mask = Parameter(torch.randn(E, device=device) * std)
    #     else:
    #         assert False