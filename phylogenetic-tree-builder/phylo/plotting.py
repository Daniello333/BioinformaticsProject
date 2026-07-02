"""Simple dendrogram rendering for a TreeNode, using matplotlib line segments."""
from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .tree import TreeNode


def _layout(node: TreeNode, depth: float, y_counter: list[int], positions: dict[int, tuple[float, float]]):
    """Assign (x=depth, y=vertical slot) to every node; returns this node's y."""
    node_id = id(node)
    if node.is_leaf():
        y = y_counter[0]
        y_counter[0] += 1
        positions[node_id] = (depth, y)
        return y

    child_ys = []
    for child, branch_len in node.children:
        child_y = _layout(child, depth + branch_len, y_counter, positions)
        child_ys.append(child_y)
    y = sum(child_ys) / len(child_ys)
    positions[node_id] = (depth, y)
    return y


def plot_dendrogram(root: TreeNode, out_path: str, title: str = "Phylogenetic Tree") -> None:
    positions: dict[int, tuple[float, float]] = {}
    _layout(root, 0.0, [0], positions)

    fig, ax = plt.subplots(figsize=(8, 5))

    def draw(node: TreeNode):
        x0, y0 = positions[id(node)]
        for child, _ in node.children:
            x1, y1 = positions[id(child)]
            # horizontal segment for the branch, vertical segment connecting siblings
            ax.plot([x0, x1], [y1, y1], color="black", linewidth=1)
            ax.plot([x0, x0], [y0, y1], color="black", linewidth=1)
            draw(child)
        if node.is_leaf():
            ax.text(x0 + 0.01 * max(1, positions[id(root)][0]), y0, node.name, va="center", fontsize=9)

    draw(root)
    ax.set_xlabel("Distance")
    ax.set_yticks([])
    ax.set_title(title)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
