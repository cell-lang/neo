# Hacker Cup 2021

## Qualification round

## Round 1

### Problem A1

Input:

```haskell
  string(Nat -> «'F', 'X', 'O'») -- This is actually a sequence
```

Output:

```haskell
  min_switches : Nat
```

Algorithm:

```haskell

  switch_char(:left) = '0';
  switch_char(:right) = 'X';

  other_hand(:left, :right);
  other_hand(:right, :left);

  string(0, C), H = :left | H = :right, other_hand(H, O), switch_char(H, Cs) ->
    switches(H, 0, if C == Cs then 1 else 0),
    hand(H, 0, if C == Cs then O else H);

  switches(H₀, I, N), hand(H₀, I, Hi), string(I + 1, C), switch_char(Hi, Cs) ->
    switches(H₀, I + 1, N + if C == Cs then 1 else 0),
    hand(H₀, 0, if C == Cs then O else H);

  last_idx = max(I) : string(I, _);
  min_switches = min(N) : switches(_, last_idx, N);

```

### Problem A2

Input:

```haskell
  string(Nat -> «'F', 'X', 'O'») -- This is actually a sequence
```

Output:

```haskell
  min_switches_sum : Nat
```

Algorithm:

```haskell

  string(I, _), string(J, C), C ≠ 'F', J < I, max J : I -> prev_non_F_char(I, J, C);

  string(I, Ci), prev_non_F_char(I, J, Cp), Ci ≠ F, Ci ≠ Cp -> delta_Gi(I, J + 1);
  
  string(I, _), delta_Gi(J, G), J ≤ I : I -> sum(G);

  Gi(_, G) -> Gs(sum(G));

```

With some syntactic sugar:

```haskell

  let string(I, Gi);

  string(J < I, C ≠ 'F'), max J : I -> prev_non_F_char(I, J, C);

  delta_Gi(I) = J + 1 : prev_non_F_char(I, J, Cp), Ci ≠ F, Ci ≠ Cp;

  Gi(I) = sum(delta_Gi(J)) : J ≤ I;

  Gs = sum(Gi);
```



<!-- 
llong solve(int N, string &S) {
  int i = 0;

  while (i < N && S[i] == 'F')
    i++;
  if (i == N)
    return 0;

  int last_non_F_char = S[i];
  int last_non_F_char_idx = i;
  i++;

  llong Gs = 0;
  llong Gi = 0;

  for ( ; i < N ; i++) {
    char ch = S[i];

    if (ch ≠ 'F') {
      if (ch ≠ last_non_F_char)
        Gi = (Gi + last_non_F_char_idx + 1) % MOD;
      last_non_F_char = ch;
      last_non_F_char_idx = i;
    }

    Gs = (Gs + Gi) % MOD;

    // cout << "i = " << i << ", Gi = " << Gi << ", Gs = " << Gs << endl;
  }

  return Gs;
}
 -->