# Basic Mathematical Specification

This document fixes the conventions used by the first synthetic Rust prototype.
It implements the graph-level decomposition from the proposal. It does not add
the optional simplicial-complex extension.

## Input signal

For one sample and one gene, each directed edge represents a splice junction
with nonnegative support `count`. The basic normalized edge signal is

\[
w_e = \frac{c_e}{\sum_j c_j}.
\]

The current prototype rejects groups whose total support is zero. It also
rejects duplicate directed edges and self-loops rather than silently deciding
how to combine them.

## Incidence matrix

Nodes and edges are sorted by their text identifiers so results do not depend
on input row order. For an edge from node \(u\) to node \(v\), its column in
\(B_1\) contains \(-1\) at \(u\), \(+1\) at \(v\), and zero elsewhere.

## Graph-level decomposition

The fitted potential solves

\[
\phi^* = \operatorname*{arg\,min}_{\phi}
\left\|F-B_1^T\phi\right\|_2^2.
\]

The two edge-space components are

\[
F_{\mathrm{grad}}=B_1^T\phi^*, \qquad
F_{\mathrm{cycle}}=F-F_{\mathrm{grad}}.
\]

The potential is not unique: adding a constant within a connected component
does not change the gradient. The implementation fixes the first sorted node of
each undirected connected component to potential zero and solves the remaining
dense system. This produces the same projected edge signal while keeping the
basic implementation dependency-free.

For nonzero signal, the reported cycle-space fraction is

\[
C=\frac{\left\|F_{\mathrm{cycle}}\right\|_2^2}
        {\left\|F\right\|_2^2}.
\]

## Acceptance checks

Every calculation checks, to a relative scale with base tolerance \(10^{-9}\),

\[
F \approx F_{\mathrm{grad}}+F_{\mathrm{cycle}},
\]

\[
B_1F_{\mathrm{cycle}} \approx 0,
\]

and

\[
F_{\mathrm{grad}}^TF_{\mathrm{cycle}} \approx 0.
\]

The command exits with an error if these checks fail.

## First two expected examples

`tests/fixtures/linear.tsv` is a three-node chain. A tree has no graph cycle
space, so every edge signal is a gradient and its expected cycle fraction is
zero.

`tests/fixtures/cycle.tsv` is the directed cycle
\(E_1 \rightarrow E_2 \rightarrow E_3 \rightarrow E_1\) with equal weights.
The signal is in \(\ker(B_1)\), so its expected gradient component is zero and
its expected cycle fraction is one.

This basic solver uses normal equations and dense Gaussian elimination. It is
appropriate for the tiny validation fixtures, but the numerical method must be
reviewed before large or ill-conditioned biological graphs are analyzed.

