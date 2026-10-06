# Periodic boundary conditions (§2.1): cell and edge indices are taken modulo N.
# Edge e = i + 1/2 is stored at index i, so cell i has edges i (right) and left(i, N) (left).

"Index of the right neighbour of cell `i` on a periodic grid of `N` cells."
@inline right(i, N) = i == N ? 1 : i + 1

"Index of the left neighbour of cell `i` on a periodic grid of `N` cells."
@inline left(i, N) = i == 1 ? N : i - 1
