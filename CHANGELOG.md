# Changelog

## Unreleased (v2)

Pydantic models for track constructors, hatch-vcs packaging, and pixi-based CI.

### Breaking

- **Python 3.10+** is required (`requires-python = ">=3.10"`).
- Unknown track keyword arguments now raise `pydantic.ValidationError`
  (`extra="forbid"`). Typos that were previously ignored will fail.
- Annotation / dual-axis **layout** height is `layout_height()`, not `.height`.
  `.height` is the per-lane unit used by GenomeViewer via `layout_height()`.
- `BigWigTrack(path, "bar")` and `OverlayingTrack(..., palette)` no longer
  accept those extra **positional** arguments; pass `plot_type=` / `palette=`.
- Constructor aliases (`inward_ticks`, `transformation`,
  `draw_y_independently`) work at init time only, not on later attribute
  assignment.
- Invalid colors, filters, and `stat_method` values raise `ValidationError`
  instead of being ignored.
- File handles such as `.bam` / `.bw` are private (`_bam`, `_bw`).

### Fixed

- Uncompressed `.gtf` / `.bedpe` files are not opened with Tabix when pysam
  is installed. GtfTrack has a pandas fallback parser.
- `LogoTrack.values` keeps the A/C/G/T (or amino-acid) DataFrame when given
  a numpy array.
- `add_highlight_region` records default color/alpha so highlights draw.
- Coverage tracks apply `data_transform` and `scale`.
- GtfTrack no longer overwrites `.height` during layout.
- Docs install command is `pip install GenomeViewer`.
- PyPI publish requires a git tag whose hatch-vcs version matches the tag.
