# Reviewed synthetic expectations

These expectations were reviewed on October 12, 2026. The values are computed
from the fixture counts and checked independently of the Rust output.

| Fixture | Expected result |
|---|---|
| `linear.tsv` | cycle fraction `0` |
| `branching_tree.tsv` | cycle fraction `0` because the graph is a tree |
| `cycle.tsv` | cycle fraction `1` for equal circulation |
| `exon_skip.tsv` | cycle fraction `1/9` |
| `noise.tsv` | cycle fraction `5618/45027` (`0.1247695827`), using fixed perturbations `+3, 0, -3` from counts of `100, 100, 100` and seed note `20261012` |
| `weight_sweep.tsv`, `0.01` | cycle fraction `4489/6667` (`0.6733163342`), denominator `201` |
| `weight_sweep.tsv`, `0.05` | cycle fraction `1681/2403` (`0.6995422389`), denominator `205` |
| `weight_sweep.tsv`, `0.10` | cycle fraction `49/67` (`0.7313432836`), denominator `210` |
| `weight_sweep.tsv`, `0.20` | cycle fraction `121/153` (`0.7908496732`), denominator `220` |
| `weight_sweep.tsv`, `0.40` | cycle fraction `8/9` (`0.8888888889`), denominator `240` |
| `disconnected.tsv` | accepted; two components and cycle fraction `0` |
| `zero_support.tsv` | rejected because total support is zero |
| `duplicate_edge.tsv` | rejected because the directed edge is duplicated |
| `self_loop.tsv` | rejected because self-loops are unsupported |

The noise fixture is a fixed, reproducible perturbation of the exon-skipping
shape. The weight sweep keeps the two forward edges at `100` and changes only
the `E3 -> E1` back-splice edge to `1`, `5`, `10`, `20`, and `40`.
