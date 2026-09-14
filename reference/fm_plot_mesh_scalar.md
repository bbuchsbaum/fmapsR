# Plot a Scalar Function on a Mesh

Plot a Scalar Function on a Mesh

## Usage

``` r
fm_plot_mesh_scalar(
  x,
  values,
  faces = NULL,
  engine = c("auto", "plotly", "base"),
  title = NULL,
  palette = "Viridis",
  show_edges = FALSE,
  view = c("isometric", "xy", "xz", "yz"),
  value_range = NULL,
  shading = TRUE
)
```

## Arguments

- x:

  Mesh domain, mesh list with \`vertices\`/\`faces\`, or vertex matrix.

- values:

  Numeric vector with one scalar per vertex.

- faces:

  Optional face matrix when \`x\` is a plain vertex matrix.

- engine:

  Plot backend: \`"auto"\`, \`"plotly"\`, or \`"base"\`.

- title:

  Optional plot title.

- palette:

  Optional palette name, vector, or function.

- show_edges:

  Whether to draw triangle edges in the base backend.

- view:

  Projection used by the base backend.

- value_range:

  Optional two-value numeric range used to keep color scaling consistent
  across multiple panels.

- shading:

  Whether the base backend should apply simple directional shading to
  improve surface legibility.

## Value

Plot object. \`plotly\` backend returns a plotly widget; \`base\`
backend returns plotting metadata invisibly.
