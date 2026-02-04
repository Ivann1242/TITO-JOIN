from __future__ import annotations

from dataclasses import dataclass
from typing import List, Tuple, Set, Optional, Union


@dataclass(frozen=True)
class PairRep:
    """
    JOIN pair representation (may have star)
    - a: left side (group node id, 0..n-1)
    - b: right side (imaginary index, e.g. x, x+n, x+2n, ...)
    - star: True means (a,b)*
    """
    a: int
    b: int
    star: bool = False

    def __str__(self) -> str:
        return f"({self.a},{self.b})" + ("*" if self.star else "")


class InversionSet:
    """
    Inversion Set data structure
    
    There are n nodes, representing groups from 0 to n-1:
    - node0 can represent 0, n, 2n, ...
    - node1 can represent 1, n+1, 2n+1, ...
    - and so on
    
    Fields:
        inversions: [(a,b), (c,d), ...]  - inversion pairs
        n: int - number of groups
    
    Class Functions:
        - convert_matrix(): convert to n*n adjacency matrix, elements are lists of reachable imaginary indices
        - JOIN_matrix(): convert to JOIN_graph, preserve all imaginary indices, mark cycles with '*'
        - compute_JOIN(): list all reachable pairs, handle star special cases
    """

    def __init__(self, inversions: List[Tuple[int, int]], n: int):
        if n <= 0:
            raise ValueError("n must be positive")
        self.inversions = inversions
        self.n = n

        # caches
        self._adj_list: Optional[List[List[int]]] = None              # group-level adjacency list
        self._adj_mat: Optional[List[List[List[int]]]] = None         # n*n matrix, elements are lists
        self._join_mat: Optional[List[List[Union[List[int], str]]]] = None  # JOIN matrix (for display)
        self._join_data: Optional[List[List[List[int]]]] = None       # JOIN data (preserves all imaginary indices)
        self._star_nodes: Optional[Set[int]] = None                   # nodes with self-loops

    # ---------------- Function 1: convert_matrix ----------------
    def convert_matrix(self) -> List[List[List[int]]]:
        """
        Convert inversion_set to n*n adjacency graph/matrix
        
        Conversion rules:
        - For each pair (a, b) in inversions, left to right means reachable
        - Matrix elements are arrays/lists, recording all reachable imaginary indices
        
        Example when n=2:
        - (0,1): matrix[0][1] records 1
        - (0,3): matrix[0][1] also records 3 (since 3 % 2 = 1)
        - Final matrix[0][1] = [1, 3]
        """
        n = self.n
        # Adjacency matrix, each element is a list recording all original imaginary indices
        adj_mat: List[List[List[int]]] = [[[] for _ in range(n)] for _ in range(n)]
        
        for a, b in self.inversions:
            u = a % n  # source group
            v = b % n  # target group
            # Record original imaginary index b
            if b not in adj_mat[u][v]:
                adj_mat[u][v].append(b)
        
        self._adj_mat = adj_mat
        
        # Also build group-level adjacency list for graph algorithms
        adj_list: List[List[int]] = [[] for _ in range(n)]
        for u in range(n):
            for v in range(n):
                if adj_mat[u][v] and v not in adj_list[u]:
                    adj_list[u].append(v)
        
        self._adj_list = adj_list

        # Clear downstream caches
        self._join_mat = None
        self._reach = None
        self._star_nodes = None

        return adj_mat

    def _ensure_graph(self) -> None:
        """Ensure adjacency matrix and adjacency list are built"""
        if self._adj_list is None or self._adj_mat is None:
            self.convert_matrix()

    # ---------------- Function 2: JOIN_matrix ----------------
    def JOIN_matrix(self) -> List[List[Union[List[int], str]]]:
        """
        Convert adjacency matrix to JOIN_graph (transitive closure)
        
        Conversion rules:
        - Preserve all reachable imaginary indices (no information loss)
        - If a cycle exists (node can reach itself via a path), the node has a self-loop
        - Mark '*' at diagonal position (i,i) to indicate self-loop
        
        Example n=2, inversions=[(0,1), (0,3), (1,2)]:
        - convert_matrix: [[[], [1,3]], [[2], []]]
        - JOIN_matrix: preserve all imaginary indices after transitive closure
        - Mark '*' on diagonal when cycle exists
        """
        self._ensure_graph()
        n = self.n

        # Initialize: copy adjacency matrix
        join_data: List[List[List[int]]] = [[list(cell) for cell in row] for row in self._adj_mat]

        # Floyd-Warshall transitive closure, preserve all imaginary indices
        changed = True
        while changed:
            changed = False
            for k in range(n):
                for i in range(n):
                    for j in range(n):
                        # If i->k reachable and k->j reachable
                        if join_data[i][k] and join_data[k][j]:
                            # Add k->j imaginary indices to i->j
                            for idx in join_data[k][j]:
                                if idx not in join_data[i][j]:
                                    join_data[i][j].append(idx)
                                    changed = True

        # Identify self-loop nodes
        self._star_nodes = {i for i in range(n) if join_data[i][i]}
        self._join_data = join_data

        # Build display matrix: use '*' for diagonal self-loops
        display_mat: List[List[Union[List[int], str]]] = []
        for i in range(n):
            row: List[Union[List[int], str]] = []
            for j in range(n):
                if i == j and i in self._star_nodes:
                    row.append('*')
                else:
                    row.append(join_data[i][j])
            display_mat.append(row)

        self._join_mat = display_mat
        return display_mat

    # ---------------- Function 3: compute_JOIN ----------------
    def compute_JOIN(self) -> Set[PairRep]:
        """
        List all reachable pairs from JOIN_graph
        
        Special case handling:
        1. If the right side (target) node x has matrix (x,x) = '*' (self-loop),
           record as star representation.
           (a,b)* represents the entire congruent class: (a,b), (a,b+n), (a,b+2n)...
           So only keep the smallest imaginary index in the same group
           
        2. All nodes x with self-loops are additionally recorded as (x, x+n)*
        
        Returns:
            JOIN - set of all reachable pairs after special case processing
        """
        if self._join_data is None or self._star_nodes is None:
            self.JOIN_matrix()

        assert self._join_data is not None
        assert self._star_nodes is not None

        join_data = self._join_data
        star_nodes = self._star_nodes
        n = self.n

        join: Set[PairRep] = set()

        for a in range(n):
            for b_group in range(n):
                if not join_data[a][b_group]:
                    continue

                # Skip (a, a) direct output, handled by special case 2
                if b_group == a:
                    continue

                if b_group in star_nodes:
                    # Target group has star: (a,b)* covers entire congruent class
                    # Only keep the smallest imaginary index
                    b_min = min(join_data[a][b_group])
                    join.add(PairRep(a, b_min, star=True))
                else:
                    # Target group has no star: keep all imaginary indices
                    for b in join_data[a][b_group]:
                        join.add(PairRep(a, b, star=False))

        # Special case 2: all nodes x with self-loops are recorded as (x, x+n)*
        for x in star_nodes:
            join.add(PairRep(x, x + n, star=True))

        return join

# ------------------ Demo ------------------
if __name__ == "__main__":
    # Example 1: no cycle
    print("n=2, inversions=[(0,1), (0,3)]")
    inv = InversionSet(inversions=[(0, 1), (0, 3)], n=2)
    
    print("\nconvert_matrix:")
    for row in inv.convert_matrix():
        print(row)
    
    print("\nJOIN_matrix:")
    for row in inv.JOIN_matrix():
        print(row)
    
    print("\ncompute_JOIN:")
    for p in sorted(inv.compute_JOIN(), key=lambda x: (x.a, x.b)):
        print(p)

    # Example 2: with cycle
    print("\n" + "=" * 30)
    print("n=2, inversions=[(0,1), (0,3), (1,2)]")
    inv2 = InversionSet(inversions=[(0, 1), (0, 3), (1, 2)], n=2)
    
    print("\nconvert_matrix:")
    for row in inv2.convert_matrix():
        print(row)
    
    print("\nJOIN_matrix:")
    for row in inv2.JOIN_matrix():
        print(row)
    
    print("\ncompute_JOIN:")
    for p in sorted(inv2.compute_JOIN(), key=lambda x: (x.a, x.b)):
        print(p)
    
    # Example 3
    print("n=3, inversions=[(0,1), (1,2), (0,3)]")
    inv = InversionSet(inversions=[(0, 1),(1, 2) ,(0, 3)], n=3)
    
    print("\nconvert_matrix:")
    for row in inv.convert_matrix():
        print(row)
    
    print("\nJOIN_matrix:")
    for row in inv.JOIN_matrix():
        print(row)
    
    print("\ncompute_JOIN:")
    for p in sorted(inv.compute_JOIN(), key=lambda x: (x.a, x.b)):
        print(p)
