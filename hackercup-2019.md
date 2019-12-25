# Facebook hacker cup 2019


## Qualification round


### Problem A

```haskell
  -- input: "A..BB..B"
  lilipad(1, alpha);
  lilipad(2, empty);
  lilipad(3, empty);
  lilipad(4, beta);
  lilipad(5, beta);
  lilipad(6, empty);
  lilipad(7, empty);
  lilipad(8, beta);


  betas = |lilipad(_, beta)|;
  lilipads = |lilipad|;

  can_jump = betas < lilipads - 1 and 2 * betas < len - 1;
```


### Problem B

```haskell
  betas = |lilipad(_, beta)|;
  lilipads = |lilipad|;

  can_jump = betas < lilipads - 1 and (betas >= 2 or 2 * betas < len - 1);
```


### Problem C

```haskell

```


### Problem D

```haskell

```
