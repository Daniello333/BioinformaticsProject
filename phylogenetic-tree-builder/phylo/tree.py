"""Tree construction: UPGMA (distance + molecular clock assumption) and
Neighbor-Joining (distance, no clock assumption)."""
from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class TreeNode:
    name: str | None = None
    children: list[tuple["TreeNode", float]] = field(default_factory=list)

    def is_leaf(self) -> bool:
        return not self.children

    def to_newick(self) -> str:
        if self.is_leaf():
            return self.name or ""
        parts = [f"{child.to_newick()}:{branch_len:.6f}" for child, branch_len in self.children]
        inner = ",".join(parts)
        return f"({inner})"

    def newick(self) -> str:
        return self.to_newick() + ";"

    def leaf_names(self) -> list[str]:
        if self.is_leaf():
            return [self.name] if self.name else []
        names = []
        for child, _ in self.children:
            names.extend(child.leaf_names())
        return names


def _copy_matrix(labels: list[str], matrix: list[list[float]]) -> dict[str, dict[str, float]]:
    dist: dict[str, dict[str, float]] = {a: {} for a in labels}
    for i, a in enumerate(labels):
        for j, b in enumerate(labels):
            dist[a][b] = matrix[i][j]
    return dist


def upgma(labels: list[str], matrix: list[list[float]]) -> TreeNode:
    """Build a rooted tree using UPGMA (unweighted pair group method with arithmetic mean)."""
    dist = _copy_matrix(labels, matrix)
    nodes: dict[str, TreeNode] = {name: TreeNode(name=name) for name in labels}
    sizes: dict[str, int] = {name: 1 for name in labels}
    heights: dict[str, float] = {name: 0.0 for name in labels}
    active = list(labels)

    while len(active) > 1:
        best_pair = None
        best_dist = float("inf")
        for i in range(len(active)):
            for j in range(i + 1, len(active)):
                a, b = active[i], active[j]
                if dist[a][b] < best_dist:
                    best_dist = dist[a][b]
                    best_pair = (a, b)

        a, b = best_pair
        new_height = best_dist / 2.0
        new_name = f"({a}+{b})"
        new_node = TreeNode(
            children=[
                (nodes[a], new_height - heights[a]),
                (nodes[b], new_height - heights[b]),
            ]
        )

        new_dist: dict[str, float] = {}
        for c in active:
            if c in (a, b):
                continue
            weighted = (sizes[a] * dist[a][c] + sizes[b] * dist[b][c]) / (sizes[a] + sizes[b])
            new_dist[c] = weighted

        active = [c for c in active if c not in (a, b)]
        for c in active:
            dist[c][new_name] = new_dist[c]
            dist[new_name] = dist.get(new_name, {})
            dist[new_name][c] = new_dist[c]
        dist[new_name] = dist.get(new_name, {})

        nodes[new_name] = new_node
        sizes[new_name] = sizes[a] + sizes[b]
        heights[new_name] = new_height
        active.append(new_name)

    return nodes[active[0]]


def neighbor_joining(labels: list[str], matrix: list[list[float]]) -> TreeNode:
    """Build an unrooted tree (represented rooted at the final join) using Neighbor-Joining."""
    dist = _copy_matrix(labels, matrix)
    nodes: dict[str, TreeNode] = {name: TreeNode(name=name) for name in labels}
    active = list(labels)

    while len(active) > 2:
        n = len(active)
        r = {a: sum(dist[a][b] for b in active if b != a) for a in active}

        best_pair = None
        best_q = float("inf")
        for i in range(n):
            for j in range(i + 1, n):
                a, b = active[i], active[j]
                q = (n - 2) * dist[a][b] - r[a] - r[b]
                if q < best_q:
                    best_q = q
                    best_pair = (a, b)

        a, b = best_pair
        branch_a = 0.5 * dist[a][b] + (r[a] - r[b]) / (2 * (n - 2))
        branch_b = dist[a][b] - branch_a
        branch_a, branch_b = max(branch_a, 0.0), max(branch_b, 0.0)

        new_name = f"({a}+{b})"
        new_node = TreeNode(children=[(nodes[a], branch_a), (nodes[b], branch_b)])

        new_dist: dict[str, float] = {}
        for c in active:
            if c in (a, b):
                continue
            new_dist[c] = 0.5 * (dist[a][c] + dist[b][c] - dist[a][b])

        active = [c for c in active if c not in (a, b)]
        for c in active:
            dist[c][new_name] = new_dist[c]
            dist.setdefault(new_name, {})[c] = new_dist[c]

        nodes[new_name] = new_node
        active.append(new_name)

    # Join the final two nodes/clusters with the remaining branch length split between them.
    a, b = active
    d = dist[a][b]
    root = TreeNode(children=[(nodes[a], d / 2), (nodes[b], d / 2)])
    return root
