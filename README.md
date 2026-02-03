# TITO JOIN — README

This project computes the JOIN closure induced by an inversion set over `n` residue groups. Each group node `i` represents all indices congruent to `i (mod n)`:
`i, i+n, i+2n, ...`. Each inversion `(a, b)` induces a directed edge from `a mod n` to `b mod n`.

The output includes:
1) A `JOIN_matrix` (size `n x n`) with entries in `{0, 1, *}`:
   - `1` if `j` is reachable from `i` by a directed path of length ≥ 1
   - `*` on the diagonal `(i, i)` if there exists a directed cycle returning to `i` (a non-zero-length path from `i` back to `i`)
   - `0` otherwise

2) A JOIN set of reachable pairs with a star representation and a canonical "imaginary index lifting" rule that avoids printing low representatives that are not strictly greater than the left endpoint.

## Classes

- `PairRep(a, b, star)`: a hashable representation of an output pair, printed as `(a,b)` or `(a,b)*`.
- `InversionSet(inversions, n)`: main class implementing:
  - `convert_matrix()`
  - `JOIN_matrix()`
  - `compute_JOIN()`

## Core Definitions

Given `n` groups indexed `0..n-1` and a list of inversions `[(a,b), ...]`:

- Group-level edge:
  - `u = a mod n`, `v = b mod n`, add directed edge `u -> v`

- Reachability (path length ≥ 1):
  - `reach[i][j] = True` iff there exists a directed path from `i` to `j` with at least one edge.

- Star nodes:
  - Node `x` is a star node iff `reach[x][x] = True` (there is a non-zero-length cycle returning to `x`).

- `JOIN_matrix`:
  - `JOIN_matrix[i][j] = "1"` if `reach[i][j] = True` and `i != j`
  - `JOIN_matrix[i][i] = "*"` if `reach[i][i] = True`
  - `JOIN_matrix[i][j] = "0"` otherwise

## Output JOIN Set and Canonical Lifting

After computing `reach`, we enumerate all reachable group pairs `(a, b_group)` with length ≥ 1 and `b_group != a`.

We then "lift" the right endpoint to a canonical imaginary index `b` in the same residue class such that `b > a`:
- `b = b_group + k*n` where `k` is the smallest integer `>= 0` making `b > a`.

This guarantees:
- We never output a right endpoint `b` that is less than or equal to `a`.
- In particular, for `n=2`, if `1` can reach group `0`, we output `(1,2)` rather than `(1,0)`.

Star marking rule:
- If `b_group` is a star node, output pair is marked as starred: `(a, b)*`.
- Otherwise it is unstarred: `(a, b)`.

Special rule for cycles:
- We do not output `(x,x)*` directly.
- For every star node `x`, we add `(x, x+n)*`.

This matches the intended format (e.g., for `n=2` and edges forming a cycle, outputs include `(0,2)*` and `(1,3)*` rather than `(0,0)*` and `(1,1)*`).

## Algorithm Overview

1) Build group-level directed graph:
   - Reduce each inversion `(a,b)` to edge `(a mod n) -> (b mod n)`.

2) Compute reachability matrix `reach` (transitive closure for paths of length ≥ 1):
   - For each start node `i`, run DFS/BFS to mark all reachable nodes `j`.
   - Record `reach[i][j] = True` for all visited `j`.

3) Construct `JOIN_matrix`:
   - Fill with `"1"` wherever `reach[i][j]` is true.
   - Replace diagonal with `"*"` where `reach[i][i]` is true.

4) Construct JOIN set:
   - For each `a` and each `b_group` with `reach[a][b_group] = True` and `b_group != a`:
     - Lift `b_group` to the smallest `b > a` in that residue class.
     - Add `(a,b)*` if `b_group` is a star node, else `(a,b)`.
   - For each star node `x`, add `(x, x+n)*`.

## Correctness Sketch

We prove that the outputs match the definitions above.

Lemma 1 (Graph construction correctness).
For every input inversion `(a,b)`, the code adds the directed edge `(a mod n) -> (b mod n)`. Therefore the constructed graph is exactly the group-level reduction of the inversion set.

Lemma 2 (Reachability correctness).
For each start node `i`, DFS/BFS visits exactly the nodes reachable from `i` by a directed path of length ≥ 1. Thus `reach[i][j]` is true iff `j` is reachable from `i` (length ≥ 1).

Lemma 3 (JOIN_matrix correctness).
By Lemma 2, `reach[i][j]` is true exactly when `j` is reachable from `i`. The code writes `"1"` for all such positions and writes `"*"` on diagonal entries where `reach[i][i]` is true. Therefore `JOIN_matrix` matches the required `{0,1,*}` definition.

Lemma 4 (Canonical lifting correctness).
For any reachable group target `b_group`, the function `lift_to_imaginary_gt_a(a, b_group)` returns the smallest index `b` in the residue class `b_group (mod n)` satisfying `b > a`. Hence each JOIN output pair is a valid representative of its target group and respects the required canonical format.

Lemma 5 (JOIN set correctness with star rule).
For any reachable pair `(a, b_group)` with `b_group != a`, the code outputs exactly one canonical lifted representative `(a,b)` and marks it starred iff `b_group` is a star node. The code also adds `(x, x+n)*` for every star node `x`, and does not output `(x,x)*` directly. Therefore the JOIN set matches the project’s output rules.

Combining Lemmas 1–5 proves overall correctness.

## Theorically Time Complexity

Let `n` be the number of group nodes, and `m` be the number of unique edges in the group graph after modulo reduction.

- Building the graph: `O(len(inversions))` time, `O(n + m)` space.
- Reachability via DFS/BFS from each node:
  - `O(n*(n + m))` time in total (each search is `O(n+m)`).
  - `O(n^2)` space to store the `reach` matrix.


## Example

Example with `n=2`, inversions `[(0,1),(1,2)]`:
- Edges: `0 -> 1`, `1 -> 0` (since `2 mod 2 = 0`)
- `JOIN_matrix`:
  - `['*', '1']`
  - `['1', '*']`
- JOIN output:
  - `(0,1)*`, `(0,2)*`, `(1,2)*`, `(1,3)*`
