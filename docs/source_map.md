# Project source map

Updated November 2, 2026.

This file records where the documents that control the project live. It keeps
the repository honest about what is available locally and what still lives in
Overleaf or in the presentation files.

| Item | Location | Purpose | Status |
|---|---|---|---|
| Scientific proposal | [`assets/hodgetheorynew.pdf`](../assets/hodgetheorynew.pdf) | Source of the research question, Hodge approach, and planned deliverables | Available locally |
| Roadmap section | [`roadmap.tex`](../roadmap.tex) | Include-only LaTeX section intended for `\\input{roadmap}` | Available locally |
| Main proposal `.tex` | Overleaf project | Parent document, preamble, title page, and section ordering | Not present in this checkout |
| Execution roadmap | [`ROADMAP.md`](../ROADMAP.md) | Detailed phases, dependencies, and acceptance gates | Available locally |
| Task guide | [`delegations_and_tasks/tasks.md`](../delegations_and_tasks/tasks.md) | Beginner-readable task breakdown and status board | Available locally |
| Math specification | [`docs/math_spec.md`](math_spec.md) | Equations, orientation, numerical checks, and limits | Available locally |
| Data methods | [`docs/data_methods.md`](data_methods.md) | Reference genome, annotation, and junction conventions | Available locally |

The repository does not silently recreate the missing parent proposal. When the
Overleaf source is available, record its filename, document class, packages,
and `\\input` or `\\include` structure here before editing it. The local
`roadmap.tex` inherits the parent document's style; it is not a replacement for
the full proposal.

## LaTeX handoff checklist

Before changing the proposal, copy the parent `.tex` source into a private
working checkout or export it from Overleaf, then record:

1. the main filename and folder;
2. the document class and page geometry;
3. the packages already in use;
4. the section file that includes `roadmap.tex`;
5. the command used to compile the PDF; and
6. any remaining overfull boxes or missing references.

The current checkout can validate the include-only section's syntax only after
the parent document is supplied. No proposal rewrite is implied by this map.
