# Add a Map to an FM Network

Add a Map to an FM Network

## Usage

``` r
fm_network_add_map(
  network,
  i,
  j,
  map,
  weight = 1,
  add_reverse = !network$directed
)
```

## Arguments

- network:

  \`fm_network\` object.

- i:

  Source domain name.

- j:

  Target domain name.

- map:

  \`fm_fit\` object or matrix.

- weight:

  Edge weight.

- add_reverse:

  For undirected networks, add reverse map automatically.

## Value

Updated \`fm_network\` object.
