# Extract 3-Cycles from Network

Extract 3-Cycles from Network

## Usage

``` r
fm_network_cycles(network)
```

## Arguments

- network:

  \`fm_network\` object.

## Value

List of 3-cycles. For directed networks each cycle is ordered
\`c(i,j,k)\` satisfying edges \`i-\>j\`, \`j-\>k\`, and \`k-\>i\`.
