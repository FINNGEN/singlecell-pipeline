# Single-cell data processing pipeline in FinnGen

Single-nucleus multiome (snRNA-seq + snATAC-seq) processing and QTL discovery pipeline used in:

Kanai, M. et al. [Population-scale immune multiome atlas reveals regulatory disease mechanisms](https://doi.org/10.1038/s41586-026-11078-2). _Nature_ **659**, 729–743 (2026)

Figure generation code can be found at <https://github.com/mkanai/finngen-multiome-flagship>.

## Repository layout

- `wdl/` — WDL workflows for each pipeline.
- `docker/` — Dockerfiles for the container images referenced by the WDLs.
- `python/` — helper scripts baked into the Docker images.
- `R/` — helper scripts baked into the Docker images.

## License

MIT License

## Contact

Masahiro Kanai (<mkanai@broadinstitute.org>)
