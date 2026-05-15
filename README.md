# Single-cell data processing pipeline in FinnGen

Single-nucleus multiome (snRNA-seq + snATAC-seq) processing and QTL discovery pipeline used in:

Kanai, M. et al. [Population-scale multiome immune cell atlas reveals complex disease drivers](https://doi.org/10.1101/2025.11.25.25340489). medRxiv (2025)

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
