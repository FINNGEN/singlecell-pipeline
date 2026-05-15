# WDL workflows

Each subdirectory contains a WDL pipeline plus its input-JSON templates. `*.tasks.wdl` files define shared tasks imported by the workflow WDLs in the same directory.

## Preprocessing

- [`demultiplexing/`](demultiplexing/): pool-level genetic demultiplexing.
  - [`cumulus_demultiplexing.wdl`](demultiplexing/cumulus_demultiplexing.wdl): demultiplex multiplexed pools with Cumulus (Demuxlet/Souporcell wrappers).
- [`generate_vcf/`](generate_vcf/): donor genotype VCFs used downstream for demultiplexing and QTL mapping.
  - [`generate_vcf.wdl`](generate_vcf/generate_vcf.wdl): build a multi-sample VCF from imputed genotypes.
  - [`subset_multiplex_vcf.wdl`](generate_vcf/subset_multiplex_vcf.wdl): subset the multi-sample VCF to the donors present in each pool.
- [`salsa/`](salsa/): SALSA-based per-modality processing and fragment splitting.
  - [`salsa.rna.wdl`](salsa/salsa.rna.wdl): per-pool snRNA-seq processing (filtering, doublet calling, count matrices).
  - [`salsa.atac.wdl`](salsa/salsa.atac.wdl): per-pool snATAC-seq processing (peak calling, fragment QC).
  - [`create_atac_objects.wdl`](salsa/create_atac_objects.wdl): build ArchR/Signac objects per pool.
- [`preprocess_h5ad/`](preprocess_h5ad/): assemble and downsample integrated H5AD objects.
  - [`preprocess_h5ad.gex.wdl`](preprocess_h5ad/preprocess_h5ad.gex.wdl) / [`preprocess_h5ad.atac.wdl`](preprocess_h5ad/preprocess_h5ad.atac.wdl): per-cell-type integration into H5AD.
  - [`downsample_h5ad.wdl`](preprocess_h5ad/downsample_h5ad.wdl) / [`downsample_h5ad.CD4_CD8_T.wdl`](preprocess_h5ad/downsample_h5ad.CD4_CD8_T.wdl) / [`downsample_h5ad.CD4_CD8_T.atac.wdl`](preprocess_h5ad/downsample_h5ad.CD4_CD8_T.atac.wdl): downsample cells per donor for sensitivity runs.
- [`multilinker/`](multilinker/): peak-gene linking with `multilinker`.
  - [`multilinker.wdl`](multilinker/multilinker.wdl): compute peak-gene links per cell type.
  - [`multilinker.optimize_n_lsi.wdl`](multilinker/multilinker.optimize_n_lsi.wdl): grid-search the number of LSI components.
  - [`multilinker.downsample.wdl`](multilinker/multilinker.downsample.wdl): linker recomputation on downsampled cells.
- [`umap/`](umap/): embeddings and batch metrics.
  - [`scanpy_umap.wdl`](umap/scanpy_umap.wdl): Scanpy UMAP on the snRNA side.
  - [`snapatac2.umap.wdl`](umap/snapatac2.umap.wdl): SnapATAC2 UMAP on the snATAC side.
  - [`compute_lisi.wdl`](umap/compute_lisi.wdl): LISI batch-mixing metric.

## QTL discovery

- [`tensorqtl/`](tensorqtl/): *cis* eQTL and caQTL mapping with TensorQTL.
  - [`tensorqtl.gex.wdl`](tensorqtl/tensorqtl.gex.wdl) / [`tensorqtl.atac.wdl`](tensorqtl/tensorqtl.atac.wdl): main *cis*-QTL scans.
  - [`tensorqtl.cis.optimize_n_peer.wdl`](tensorqtl/tensorqtl.cis.optimize_n_peer.wdl) / [`...optimize_n_expression_pcs.wdl`](tensorqtl/tensorqtl.cis.optimize_n_expression_pcs.wdl) / [`...optimize_n_loco_expression_pcs.wdl`](tensorqtl/tensorqtl.cis.optimize_n_loco_expression_pcs.wdl): sweep the number of covariate factors (PEER, expression PCs, LOCO PCs).
  - [`tensorqtl.cis.optimize_n_peer.downsample.wdl`](tensorqtl/tensorqtl.cis.optimize_n_peer.downsample.wdl) / [`tensorqtl.cis.optimize_n_expression_pcs.downsample.wdl`](tensorqtl/tensorqtl.cis.optimize_n_expression_pcs.downsample.wdl) / [`tensorqtl.downsample.gex.wdl`](tensorqtl/tensorqtl.downsample.gex.wdl): downsampled version for sensitivity analyses.
  - [`tensorqtl.export_zfiles.wdl`](tensorqtl/tensorqtl.export_zfiles.wdl): export per-locus z-files for downstream fine-mapping/SuSiE.
- [`saige_qtl/`](saige_qtl/): single-cell QTL mapping with SAIGE-QTL.
  - [`saige_qtl.step1.wdl`](saige_qtl/saige_qtl.step1.wdl) / [`saige_qtl.step1.downsample.wdl`](saige_qtl/saige_qtl.step1.downsample.wdl): null model fitting.
  - [`saige_qtl.step2.wdl`](saige_qtl/saige_qtl.step2.wdl): single-variant association.
  - [`saige_qtl.export_zfiles.wdl`](saige_qtl/saige_qtl.export_zfiles.wdl): export per-locus z-files.
- [`mashr/`](mashr/): cross-cell-type effect sharing.
  - [`mashr.wdl`](mashr/mashr.wdl): MASHR on QTL summary statistics across cell types.

## Fine-mapping, colocalization & integration

- [`susie_ldcov/`](susie_ldcov/): covariate-adjusted SuSiE fine-mapping.
  - [`susie_ldcov.wdl`](susie_ldcov/susie_ldcov.wdl): per-locus SuSiE with covariate-adjusted LD (`ldcov`).
- [`smr/`](smr/): Summary-data-based Mendelian Randomization.
  - [`munge_tensorqtl.wdl`](smr/munge_tensorqtl.wdl) / [`munge_saige_qtl.wdl`](smr/munge_saige_qtl.wdl) / [`munge_fg_sumstats.wdl`](smr/munge_fg_sumstats.wdl): prepare exposure/outcome summary statistics.
  - [`smr.wdl`](smr/smr.wdl): run SMR.

## Heritability

- [`ldsc/`](ldsc/): stratified LD-score regression.
  - [`compute_ldscore.wdl`](ldsc/compute_ldscore.wdl): compute per-annotation LD scores.
  - [`sldsc.wdl`](ldsc/sldsc.wdl): run stratified LDSC against GWAS summary statistics.
- [`mesc/`](mesc/): mediated expression-score regression.
  - [`compute_expscore.indiv.wdl`](mesc/compute_expscore.indiv.wdl): individual-level expression-score computation.
  - [`mesc.wdl`](mesc/mesc.wdl) / [`mesc.bulk.wdl`](mesc/mesc.bulk.wdl): run MESC (single-cell / bulk).
  - [`munge_mesc.tensorqtl.wdl`](mesc/munge_mesc.tensorqtl.wdl) / [`munge_mesc.saige_qtl.wdl`](mesc/munge_mesc.saige_qtl.wdl) / [`munge_mesc.pqtl.wdl`](mesc/munge_mesc.pqtl.wdl) / [`munge_mesc.meta.wdl`](mesc/munge_mesc.meta.wdl): prepare TensorQTL / SAIGE-QTL / pQTL / meta summary statistics for MESC.
