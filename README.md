# Neo programming language


### Sorting, no duplicates

```haskell
  xs(X), xs(X'), X' < X : X -> ord(X, |X'|);
```


### Breadth-first search

```haskell
  distance(start, 0);
  distance(N₀, I₀), I₀ < I, edge(N₀, N₁), not (distance(N₁, I₁), I₁ < I) -> distance(N₁, I);
```

Syntactically sugared version:

```haskell
  distance(start) = 0;
  loop I {
    distance(N₀, _), edge(N₀, N₁), not distance(N₁, _) -> distance(N₁, I);
  }
```

Desugared version:

```haskell
  distance(start, 0);
  _NNNN_distance(start, 0, 0);

  _NNNN_distance(N₀, _, _I₀), _I₀ < I,
  not (_NNNN_distance(N₁, _, _I₁), _I₁ < I),
  edge(N₀, N₁) ->
    distance(N₁, I), _NNNN_distance(N₁, I, I);
```


### Depth-first search

```haskell
  order(start, 0);
  order(V₀, I₀), I₀ < N, edge(V₀, V₁), not (order(V₁, I₁), I₁ < N), max I₀, min V₁ ->
    order(V₁, N);
```

With minimal syntactic sugar:

```haskell
  order(start) = 0;
  order(V₀, I₀ < N), edge(V₀, V₁), not order(V₁, I₁ < N), max I₀, min V₁ -> order(V₁) = N;
```


With full syntactic sugar:

```haskell
  order(start) = 0;
  loop I {
    order(V₀, I₀), edge(V₀, V₁), not order(V₁, _), max I₀, min V₁ -> order(V₁) = I;
  }
```


### Heapsort (no duplicates)

```haskell
  xs(X), min X -> order(X, 0);
  xs(X₁), not (order(X₁, I₁), I₁ < N), min X₁ -> order(X₁, N);
```

With minimal syntactic sugar:

```haskell
  xs(X), min X -> order(X, 0);
  xs(X₁), not order(X₁, I₁ < N), min X₁ -> order(X₁, N);
```

With minimal syntactic sugar:

```haskell
  xs(X), min X -> order(X, 0);
  loop I {
    xs(X₁), not order(X₁, _), min X₁ -> order(X₁) = I;
  }
```

### Topological sort (no cycles)

Input:

```haskell
  vertex(V);
  edge(V₀, V₁) -> vertex(V₀), vertex(V₁);
```

Output:

```haskell
  layer(V -> L) -> vertex(V);
```

Algorithm:

```haskell
  vertex(V), not edge(_, V) -> layer(V, 0);
  layer(V₀, I₀), I₀ < N, edge(V₀, V₁), not (layer(V₁, I₁), I₁ < N) -> layer(V₁, N);
```

With minimal syntactic sugar:

```haskell
  vertex(V), not edge(_, V) -> layer(V, 0);
  layer(V₀, I₀ < N), edge(V₀, V₁), not layer(V₁, I₁ < N) -> layer(V₁, N);
```

With full syntactic sugar:

```haskell
  vertex(V), not edge(_, V) -> layer(V, 0);
  loop N {
    layer(V₀, I₀), edge(V₀, V₁), not layer(V₁, _) -> layer(V₁, N);
  }
```


### Edmonds-Karp algorithm

Input:

```haskell
  -- Input
  vertex(V);
  edge(U, V, C);
  source;
  sink;
```

State:

```haskell
  -- Edge (U, V) is part of the N-th augmenting path
  path(U, V, N);

  -- Capacity of the N-th augmenting path
  capacity(N, C);
```

Algorithm:

```haskell
  -- Adding inverse edges
  edge(U, V, C) -> edge'(V, U, C);
  edge(U, V, _), not edge(V, U, _) -> edge'(V, U, 0);

  -- Positive flow after N steps
  edge'(V₀, V₁, _), not path(V₀, V₁, I), I < N -> flow(V₀, V₁, N, 0);
  path(V₀, V₁, I), I < N, capacity(I, C) : V₀, V₁, N -> flow(V₀, V₁, N, sum(C));

  -- Residual edges and their capacity after N steps
  edge'(V₀, V₁, C), flow(V₀, V₁, N, F₀₁), flow(V₁, V₀, N, F₁₀), F = F₀₁ - F₁₀, C > F ->
    res_edge(V₀, V₁, N, C - F);

  -- Breadth-first search at Nth step

  -- If after N steps there's a residual edge between source and V, then
  -- there's also a residual path of length 1 between them
  res_edge(source, V, N, _) -> res_path(source, V, N, 1);

  -- If after N steps there's a residual path to V₀ whose length is less than L,
  -- a residual edge between V₀ and V₁ and furthermore there's no residual path
  -- to V₁ whose length is less than L, then we've found a residual path to V₁
  -- whose length is L
  res_path(_, V₀, N, L₀), L₀ < L, res_edge(V₀, V₁, N, _),
    not (res_path(_, V₁, N, L₁), L₁ < L) : V₁->
      res_path(any(V₀), V₁, N, L);

  -- If V₁ is either the sink or part of the Nth path from sink to source
  -- and V₀ is its designated predecessor in the BFS, then the egde between
  -- V₀ and V₁ is part of the Nth path from sink to source
  (path(V₁, _, N) | V₁ = sink), res_path(V₀, V₁, N, _) -> path(V₀, V₁, N);

  -- The capacity of the Nth augmenting path is the minimum capacity of its edges
  path(U, V, N), res_edge(U, V, N, C) : N -> capacity(N, min(C));

  -- Maximum network flow is the sum of the capacities of all augmenting paths
  capacity(_, C) -> max_flow(sum(C));
```

Syntactically sugared version:

```haskell
  -- Adding inverse edges
  edge(U, V, C) -> edge'(V, U, C);
  edge(U, V, _), not edge(V, U, _) -> edge'(V, U, 0);

  loop N : path, capacity {
    -- Positive flow after N steps is the sum of the capacities
    -- of all previous paths that include a given edge
    flow(V₀, V₁) += capacity(_) : path(V₀, V₁, N₀₁ < N);

    -- Residual edges and their capacity after N steps
    res_edge(V₀, V₁) = C - F : edge'(V₀, V₁, C), F = flow(V₀, V₁) - flow(V₁, V₀) if C > F;

    -- Breadth-first search at Nth step
    res_path(source, V) :- res_edge(source, V, _);
    loop I {
      res_path(_, V₀), res_edge(V₀, V₁, _), not res_path(_, V₁) : V₁ -> res_path(any(V₀), V₁)
    }

    -- If V₁ is either the sink or part of the Nth path from sink to source
    -- and V₀ is its designated predecessor in the BFS, then the egde between
    -- V₀ and V₁ is part of the Nth path from sink to source
    path(V₀, V₁, N) :- (path(V₁, _, N) | V₁ = sink), res_path(V₀, V₁);

    -- The capacity of the Nth augmenting path is the minimum capacity of its edges
    capacity(N) = min(C) : path(V₀, V₁, N), res_edge(U, V, C);
  }

  -- Maximum network flow is the sum of the capacities of all augmenting paths
  max_flow = sum(C) : capacity(_, C);
```


<br><br><br><br><br><br><br><br><hr>


### Selection sort

```haskell
  xs(X), min X -> ys(0, X);
  ys(I, Y), xs(X), X > Y, min X : I -> ys(I+1, X);
```


### Heapsort (no duplicates)

```haskell
  xs(X), not ~ys(_, X), min X -> ys(|~ys|, X);
```

Syntactically sugared versions:

```haskell
  xs(X), not ~ys(_, X), min X -> ys(~, X);
```

```haskell
  xs(X), ~X, min X -> ys(~, X);
```

```haskell
  xs(~X), min X -> ys(~, X);
```

```haskell
  xs(X), min X -> ys(~, ~X);
```


### Heapsort (stable, retains duplicates)

```haskell
  xs(I, X), not ~ixs(I), min X I -> ys(|~ys|, X), ixs(I);
```

Syntactically sugared versions:

```haskell
  xs(I, X), ~I, min X I -> ys(~, X);
```

```haskell
  xs(~I, X), min X I -> ys(~, X);
```


### Sorting by key

```haskell
  xs(~I, X), min k(X) I -> ys(~, X);
```


### Dijkstra's shortest path algorithm

Input:

```haskell
  node(N);
  edge(N₁, N₂) -> node(N₁), node(N₂);
```

Output:

```haskell
  distance(N₁, N₂, W) -> node(N₁), node(N₂);
```

Distances from a given start node `s`:

```haskell
  distance(s, 0);
  distance(N₀, W₀), edge(N₀, N₁, W₁), not ~distance(N₁, _), W = W₀ + W₁, min W -> distance(N₁, W);
```

With syntactic sugar:

```haskell
  distance(s) = 0;
  distance(N₁) <- W : distance(N₀, W₀), edge(N₀, N₁, W₁), W = W₀ + W₁, min W;
```

```haskell
  distance(s) = 0;
  distance(N₁) := W : distance(N₀, W₀), edge(N₀, N₁, W₁), W = W₀ + W₁, min W;
```

All distances between all pairs of nodes:

```haskell
  node(S) -> distance(S, S, 0);
  distance(S, N₀, W₀), edge(N₀, N₁, W₁), not ~distance(S, N₁, _), W = W₀ + W₁, min W : S ->
    distance(S, N₁, W);
```

With syntactic sugar:

```haskell
  node(S) -> distance(S, S, 0);
  distance(S, N₀, W₀), edge(N₀, N₁, W₁), S ~ N₁, W = W₀ + W₁, min W : S ->
    distance(S, N₁, W);
```

```haskell
  distance(N, N) = 0 : node(N);

  distance(~S, N, W₀), edge(N, ~T, W₁), W = W₀ + W₁, min W : S -> distance(S, T, W);
```

```haskell
  distance(N, N) = 0 : node(N);

  distance(S, N, W₀), edge(N, T, W₁), W = W₀ + W₁, min W : S -> distance(S, T -> W);
```


### Topological sort (no cycles)

Input:

```haskell
  node(N);
  edge(N₁, N₂) -> node(N₁), node(N₂);
```

Output:

```haskell
  level(N -> L) -> node(N);
```

Algorithm:

```haskell
  node(S), not edge(S, _) -> level(S, 0);
  edge(S, T), L = level(T) : S -> level(S, 1 + max(L));
```

Syntactically sugared versions:

```haskell
  level(S) = 0 : node(S), not edge(S, _);
  level(S) = 1 + max(L) : edge(S, T), L = level(T);
```

```haskell
  level(S) = 1 + max(L) | 0 : node(S), edge(S, T), L = level(T);
```

Wrong:

```
  edge(S, T), level(T, L) : S -> level(S, 1 + max(L));
```


### Depth-first search

Starting from node `n`:

```haskell
  index(n, 0);
  index(N, I), edge(N, N'), not ~index(N', _), max I ~> index(N', |~index|);
```

For all nodes of an arbitrary graph:

```haskell
  index(N, I), edge(N, N'), not ~index(N', _), max I || node(N), not ~index(N, _) ~>
    index(N', |~index|);
```

Syntactically sugared version:

```haskell
  index(N', I), edge(N', ~N), max I || node(~N) ~> index(N, ~);
```


### Tarjan's algorithm

```haskell
  --- Depth-first search
  index(V, I), edge(V, ~W), max I || vertex(~W) ~> index(W, ~);

  --- If lower_link(V) == index(V), then V is the first discovered vertex of a cycle
  lower_link(V) = min(if index(W) <= index(V) then index(W) else lower_link(W)) : edge(V, W) | W = V;

  --- lowest_link(V) is the lowest-indexed vertex that is reachable from V
  lowest_link(V) = V : lower_link(V) == index(V);
  lowest_link(V) = W : edge(V, W), min index(lowest_link(W));

  --- If V belongs to a cycle, cycle_repr(V) is the lowest-indexed vertex of it
  cycle_repr(V, V) :- lower_link(V) == index(V);
  cycle_repr(V, R) :- lowest_link(V, R), edge(W, V), cycle_repr(W, R);

  --- Representative of the strongly connected component V belongs to
  repr(V) = if cycle_repr(V, _) then cycle_repr(V, W) else V : vertex(V);
```


### Edmonds-Karp algorithm

Input:

```haskell
  -- Input
  vertex(V);
  edge(U, V, C);
  source;
  sink;
```

State:

```haskell
  -- Edge (U, V) is part of the I-th augmenting path
  path(U, V, I);

  -- Capacity of the I-th augmenting path
  capacity(I, C);
```

Algorithm:

```haskell
  -- Adding inverse edges
  edge(U, V, C), not ~edge(V, U, _) -> edge(V, U, 0);

  -- Positive flow
  flow(U, V) += capacity(I) : path(U, V, I);

  -- Single iteration of the algorithm
  atomic flow -> path, capacity {
    -- Residual edges and their capacity
    edge(U, V, C), F = ~flow(U, V) - ~flow(V, U), C > F -> res_edge(U, V, C - F);

    -- Breadth-first search
    res_edge(source, V, _) -> res_path(source, V);
    res_path(_, V), res_edge(V, W, _), not ~res_path(_, W) : W -> res_path(any(V), W);

    -- Index of the next augmenting path
    idx = 1 + max(I) | 0 : ~path(_, _, I);

    -- Tracing the new augmenting path from the sink back to the source
    path(U, V, idx) :- res_path(U, V), (path(V, _, idx) | V == sink);

    -- Capacity of the new augmenting path
    capacity(idx) = min(C) : path(U, V, idx), res_edge(U, V, C);
  }

  -- Maximum network flow is the sum of the capacities of all augmenting paths
  max_flow += C : capacity(_, C);
```

With some syntactic sugar:

```haskell
  -- Adding inverse edges
  edge(U, V, C), not ~edge(V, U, _) -> edge(V, U, 0);

  -- Positive flow
  flow(U, V) += capacity(I) : path(U, V, I);

  -- Single iteration of the algorithm
  atomic flow -> path, capacity @ idx {
    -- Residual edges and their capacity
    edge(U, V, C), F = ~flow(U, V) - ~flow(V, U), C > F -> res_edge(U, V, C - F);

    -- Breadth-first search
    res_edge(source, V, _) -> res_path(source, V);
    res_path(_, V), res_edge(V, W, _), not ~res_path(_, W) : W -> res_path(any(V), W);

    -- Tracing the new augmenting path from the sink back to the source
    path(U, V, idx) :- res_path(U, V), (path(V, _, idx) | V == sink);

    -- Capacity of the new augmenting path
    capacity(idx) = min(C) : path(U, V, idx), res_edge(U, V, C);
  }

  -- Maximum network flow is the sum of the capacities of all augmenting paths
  max_flow += C : capacity(_, C);
```

Possible alternatives:

```haskell
  -- Adding inverse edges
  edge(U, V, C), not ~edge(V, U, _) -> edge(V, U, 0);

  -- Positive flow
  flow(U, V) += capacity(I) : path(U, V, I);

  -- Residual edges and their capacity
  edge(U, V, C), F = flow(U, V) - flow(V, U), C > F -> res_edge(U, V, C - F);

  -- Breadth-first search
  res_edge(source, V, _) -> res_path(source, V);
  res_path(_, V), res_edge(V, W, _), not ~res_path(_, W) : W -> res_path(any(V), W);

  -- Index of the next augmenting path
  idx = 1 + max(I) | 0 : path(_, _, I);

  -- Single iteration of the algorithm
  atomic res_edge, res_path, idx -> path, capacity {
    -- Tracing the new augmenting path from the sink back to the source
    path(U, V, ~idx) :- ~res_path(U, V), (path(V, _, ~idx) | V == sink);

    -- Capacity of the new augmenting path
    capacity(~idx) = min(C) : path(U, V, ~idx), ~res_edge(U, V, C);
  }

  -- Maximum network flow is the sum of the capacities of all augmenting paths
  max_flow += C : capacity(_, C);
```

```haskell
  -- Adding inverse edges
  edge(U, V, C), not ~edge(V, U, _) -> edge(V, U, 0);

  -- Positive flow
  flow(U, V) += capacity(J) : path(U, V, J);

  -- Residual edges and their capacity
  edge(U, V, C), F = flow(I, U, V) - flow(I, V, U), C > F -> res_edge(I, U, V, C - F);

  -- Breadth-first search
  res_edge(source, V, _) -> res_path(source, V);
  res_path(_, V), res_edge(V, W, _), not ~res_path(_, W) : W ~> res_path(V, W);

  -- Tracing the new augmenting path from the sink back to the source
  path'(U, V) :- res_path(U, V), (path'(V, _) | V == sink);

  -- Capacity of the new augmenting path
  capacity' = min(C) : path'(U, V), res_edge(U, V, C);

  -- Inserting the new path into the knowledge base
  -- This is the tricky part
  ~path'(U, V), I = |~capacity| -> path(U, V, I), capacity(I, capacity');

  -- Maximum network flow is the sum of the capacities of all augmenting paths
  max_flow += C : capacity(_, C);
```


Is there any way to make this one work?

```haskell
  -- Adding missing inverse edges
  edge(U, V, C), not ~edge(V, U, _) -> edge(V, U, 0);

  -- Positive flow
  flow(I, U, V) += capacity(J) : I ≥ 0, path(U, V, J), 0 ≤ J < I;

  -- Residual edges and their capacity
  edge(U, V, C), F = flow(I, U, V) - flow(I, V, U), C > F -> res_edge(I, U, V, C - F);

  -- Breadth-first search
  res_edge(I, source, V, _) -> res_path(I, source, V);
  res_path(I, _, V), res_edge(I, V, W, _), not ~res_path(I, _, W) : W -> res_path(I, any(V), W);

  -- Tracing the new augmenting path from the sink back to the source
  path(U, V, I) :- res_path(I, U, V), (path(V, _, I) | V == sink);

  -- Capacity of the new augmenting path
  capacity(I) = min(C) : path(U, V, I), res_edge(I, U, V, C);

  -- Maximum network flow is the sum of the capacities of all augmenting paths
  max_flow += C : capacity(_, C);
```








<br><br><br><br><br><br><br><br><hr>

## GARBAGE GARBAGE GARBAGE

### Clustering

```
  path(V, W) :- edge(V, W) | (path(V, V'), edge(V', W));
  close(V, W) :- path(V, W), path(W, V);

  close(V, W), min W : V -> cluster(V, W);
```

```
  path(V, W) :- edge(V, W) | (path(V, V'), edge(V', W));

  close(V, V) :- vertex(V);
  close(V, W) :- path(V, W), path(W, V);

  cluster(V) = min W : close(V, W);
```


```
  path(V, W) :- edge(V, W) || path(V, V'), edge(V', W);

  close(V, W) :- vertex(V), V = W || path(V, W), path(W, V);

  cluster(V) = min W : close(V, W);
```

```
  edge(V₀, V₁) -> path(V₀, V₁);
  path(V₀, V₁), path(V₁, V₂) -> path(V₀, V₂);
```

### Topological sort for cyclical graphs

```
  edge(N₀, N₁) -> reachable(N₀, N₁);
  reachable(N₀, N₁), edge(N₁, N₂) -> reachable(N₀, N₂);


```

```
  reachable(N₀, N₁) :- edge(N₁, N₀) | (reachable(N, N₁), edge(N, N₀));

  cycle(N₀, N₁) :- reachable(N₀, N₁), reachable(N₁, N₀);

  cluster(N) = min(N, min(N')) | N : cycle(N, N');

  edge(N₀, N₁) -> cluster_edge(cluster(N₀), cluster(N₁));

  level(S) = 1 + max(L) | 0 : cluster(_, S), cluster_edge(S, T), L = level(T);

  level(N) = level(cluster(N));
```


### Lowest reachable node

```
  edge(N₀, N₁) -> reachable(N₁, N₀);
  reachable(N₁, N₀), edge(N₁, N₂) -> reachable(N₂, N₀);

  reachable(N', N), min N' : N -> lowest(N, N);
```



### Depth-first search

```
  next_idx := 1 + max(I) | 0 : ~index(_, I);

  index(N) <~ next_idx : index(N', I), edge(N', N), max I || node(N);
```

```
  next_idx := 1 + max(I) | 0 : ~index(_, I);

  index(N) = next_idx : index(N', I), edge(N', N), max I, any N;
  index(N) = next_idx : node(N), any N;
```

Procedural version:

```
  do {
    index = 0;
    for n <- node
      visit(n)

    visit(s) {
      if ~order(s) {
        yield order(s, index);
        index = index + 1;
        for t <- edge(s, ?)
          visit(t)        
      }
    }
  }
```

### Topological sort

Alternative versions:

```
  node(S), not edge(S, _) -> level(S, 0);
  node(S), not (edge(S, T), not level(T, _)) -> ready(S);
  ready(S), edge(S, T), level(T, L) : S -> level(S, 1 + sum(L));
```

```
  for L = 0..
    node(S), not (edge(S, T), not level(T, _)) -> level(S, L);
```


### Sorting (no duplicates)

```
  xs(_, X),            min X -> ys(0,   X);
  xs(_, X), X > ys(I), min X -> ys(I+1, X);
```

```
  xs(_, X), I = 0,       min X -> ys(I, X);
  xs(_, X), X > ys(I-1), min X -> ys(I, X);
```

```
  xs(_, X), I = 0 | X > ys(I-1), min X -> ys(I, X);
```

```
  xs(_, X), min X -> ys(I, X):
    I = 0 |
    X > ys(I-1);
```

```
  let xs(_, X):
    ys(0)   = min X;
    ys(I+1) = min X : X > ys(I);    
```

```
  xs(_, X):
    ys(0)   = min X;
    ys(I+1) = min X : X > ys(I);    
```

```
  ys(I) = min X : xs(_, X),
    I = 0 |
    X > ys(I-1);
```

```
  ys(I) = min X : xs(_, X), I = 0 | X > ys(I-1);
```

```
  xs(_, X), I = 0 | X > ys(I-1), min X -> ys(I, X);
```

### Sorting (stable, retains duplicates)

Syntactically sugared versions:

```
  xs(K, X), min X K -> ys(0, X), js(0, K);
  js(I, J), xs(K, X), X K > ys(I) J, min X K -> ys(I+1, X), js(I+1, K);
```

```
  xs(K, X), min X K -> ys(I, X), js(I, K):
    I = 0;
    ys(I-1, Y), js(I-1, J), X K > Y J;
```

```
  xs(J, X), min X J -> ys(I, X), js(I, J):
    I = 0;
    X J > ys(I-1) js(I-1);
```

```
  js(0)   = K : xs(K, X), min X K;
  js(I+1) = I : xs(K, X), X K > ys(I) js(I), min X K;
  
  ys(I)    = xs(js(I));
```

```
  let xs(K, X):
    js(0)   = K : min X K;
    js(I+1) = K : X K > ys(I) js(I), min X K;
  ys(I) = xs(js(I));
```

```
  js(I) = J : min X J, xs(J, X), I = 0 | X J > ys(I-1) js(I-1);
  ys(I) = xs(js(I));
```

```
  js(I) = J : min X J, xs(J, X),
    I = 0 |
    X J > ys(I-1) js(I-1);

  ys(I) = xs(js(I));
```

```
  ys(I) = X, js(I) = J : min X J, xs(J, X), I = 0 | X J > ys(I-1) js(I-1);
```

```
  xs(J, X), I = 0 | X J > ys(I-1) js(I-1), min X J -> ys(I, X), js(I, J);
```

```
  xs(J, X), I = 0 | X J > ys(I-1) js(I-1), min X J -> ys(I) = X, js(I) = J;
```


<!-- ₀₁₂ -->
