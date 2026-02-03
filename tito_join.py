from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple, Set, Optional


@dataclass(frozen=True)
class PairRep:
    """
    JOIN pair representation(may have star)
    - a: left side (group node id, 0..n-1)
    - b: right side (imaginary index, eg. x+n, x+2n, ...)
    - star: True 表示 (a,b)*
    """
    a: int
    b: int
    star: bool = False

    def __str__(self) -> str:
        return f"({self.a},{self.b})" + ("*" if self.star else "")


class InversionSet:
    """
    inversion_set:
      fields:
        inversions: [(a,b), (c,d), ...]  
        n: int                         

    functions:
      - convert_matrix(): adjacency matrix (n*n bool), adjacency list
      - JOIN_matrix(): transitive closure matrix with '0'/'1'/'*'
            '1' if reachable with path length>=1,
            '*' on diagonal if i can reach itself with path length>=1 (cycle),
            otherwise '0'
      - compute_JOIN(): list all reachable pairs with special star rules and
            imaginary-index lifting (so results match your desired format)
    """

    def __init__(self, inversions: List[Tuple[int, int]], n: int):
        if n <= 0:
            raise ValueError("n must be positive")
        self.inversions = inversions
        self.n = n

        # caches
        self._adj_list: Optional[List[List[int]]] = None
        self._adj_mat: Optional[List[List[bool]]] = None
        self._join_mat: Optional[List[List[str]]] = None
        self._reach: Optional[List[List[bool]]] = None  # reachability (len>=1)
        self._star_nodes: Optional[Set[int]] = None     # nodes with reach[i][i]=True

    # ---------------- Function 1 ----------------
    def convert_matrix(self) -> List[List[bool]]:
        """
        convert inversion_set to group-level adjacent graph/matrix
        rule: each (a,b) generate edge (a mod n) -> (b mod n), or Matrix(a,b)=1
        """
        n = self.n
        adj_mat = [[False] * n for _ in range(n)]
        adj_list = [[] for _ in range(n)]

        for a, b in self.inversions:
            u = a % n
            v = b % n
            if not adj_mat[u][v]:
                adj_mat[u][v] = True
                adj_list[u].append(v)

        self._adj_mat = adj_mat
        self._adj_list = adj_list

        # invalidate downstream caches
        self._join_mat = None
        self._reach = None
        self._star_nodes = None

        return adj_mat

    def _ensure_graph(self) -> None:
        if self._adj_list is None or self._adj_mat is None:
            self.convert_matrix()

    # ---------------- Reachability (len >= 1) ----------------
    def _reachable_from_len_ge_1(self, s: int) -> Set[int]:
        """
        all reachable from s
        """
        self._ensure_graph()
        g = self._adj_list
        visited = [False] * self.n
        stack = [s]

        while stack:
            v = stack.pop()
            for w in g[v]:
                if not visited[w]:
                    visited[w] = True
                    stack.append(w)

        return {i for i, ok in enumerate(visited) if ok}

    # ---------------- Function 2 ----------------
    def JOIN_matrix(self) -> List[List[str]]:
        """
        output n*n JOIN_matrix:
          - reachable (len>=1) => '1'
          - diagonal reachable (cycle) => '*'
          - otherwise => '0'
        """
        self._ensure_graph()
        n = self.n

        reach = [[False] * n for _ in range(n)]
        for i in range(n):
            for j in self._reachable_from_len_ge_1(i):
                reach[i][j] = True

        join_mat = [["0"] * n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if reach[i][j]:
                    join_mat[i][j] = "1"
            if reach[i][i]:
                join_mat[i][i] = "*"

        self._reach = reach
        self._join_mat = join_mat
        self._star_nodes = {i for i in range(n) if reach[i][i]}

        return join_mat

    # ---------------- Helper: lift b_group to imaginary index > a ----------------
    def _lift_to_imaginary_gt_a(self, a: int, b_group: int) -> int:
        """
        b_group in [0..n-1]
        return min imaginary index b' in the group, which satisfy b' > a
          b' = b_group + k*n, k>=0, which b' > a
        """
        n = self.n
        if b_group > a:
            return b_group
        k = ((a + 1) - b_group + n - 1) // n  # ceil((a+1-b)/n)
        return b_group + k * n

    # ---------------- Function 3 ----------------
    def compute_JOIN(self) -> Set[PairRep]:
        """
        - based on group-level reachability(len>=1)
        - if object group b 是 star node(JOIN_matrix[b][b]=='*'), then (a, lifted_b)*
        - don't output (a,a)* only follow special case 2 
        - special case 2: for all star node x, (x, x+n)*

        """
        if self._reach is None or self._star_nodes is None:
            self.JOIN_matrix()

        assert self._reach is not None
        assert self._star_nodes is not None

        reach = self._reach
        star_nodes = self._star_nodes

        join: Set[PairRep] = set()

        # enumerate reachable pairs at group-level, then lift output
        for a in range(self.n):
            for b_group in range(self.n):
                if not reach[a][b_group]:
                    continue

                # 你期望的输出不包含 (a,a)*；自环交给 special case 2
                if b_group == a:
                    continue

                b = self._lift_to_imaginary_gt_a(a, b_group)
                join.add(PairRep(a, b, star=(b_group in star_nodes)))

        # special case 2: add (x, x+n)*
        for x in star_nodes:
            join.add(PairRep(x, x + self.n, star=True))

        return join


# ------------------ Demo ------------------
if __name__ == "__main__":
    # Your example: n=2, inversions = [(0,1),(1,2)]
    inv = InversionSet(inversions=[(0, 1), (1, 2), (2, 3)], n=3)

    print("Adjacency matrix:")
    A = inv.convert_matrix()
    for row in A:
        print(row)

    print("\nJOIN matrix (* on diagonal means in a cycle):")
    Jm = inv.JOIN_matrix()
    for row in Jm:
        print(row)

    join = inv.compute_JOIN()
    print("\nJOIN (formatted):")
    for p in sorted(join, key=lambda x: (x.a, x.b, x.star)):
        print(str(p))
