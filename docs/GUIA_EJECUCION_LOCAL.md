# K-CHOPORE v3.0 — Guia de Ejecucion Local

## Resumen del Estado

```
Tu maquina (Windows + Git Bash):
  [x] Docker 27.0.3 + imagen k-chopore:latest (22.7 GB)
  [x] Python 3.13, Snakemake 9.16, pytest 9.0
  [x] results_local/ con resultados v2: FLAIR counts, DESeq2, NanoPlot, etc.
  [x] 136 tests pasando
  [ ] Reference genome TAIR10 (esta en el servidor)
  [ ] BAM alineados (estan en el servidor, ~50 GB)
  [ ] FAST5 raw signal (en NAS, ~644 GB)
```

---

## MODO A: Ejecutar Scripts v3 con Datos Locales (SIN Docker)

Esto funciona AHORA con lo que ya tienes.

### A.1 Ejecutar lncRNA Consensus (demo con mock data)

```bash
cd /c/Users/geope/Downloads/K-CHOPORE-master\ \(2\)/K-CHOPORE-master

python scripts/lncrna_consensus.py \
    --transdecoder tests/toy_data/mock_results/lncrna/transdecoder/longest_orfs.pep \
    --feelnc-gtf tests/toy_data/mock_results/lncrna/feelnc/candidate_lncRNA.gtf \
    --feelnc-class tests/toy_data/mock_results/lncrna/feelnc/candidate_lncRNA.gtf \
    --cpc2 tests/toy_data/mock_results/lncrna/cpc2/cpc2_results.txt \
    --cpat tests/toy_data/mock_results/lncrna/cpat/cpat_results.tsv \
    --candidates-gtf tests/toy_data/mock_results/lncrna/candidates_merged.gtf \
    --counts-matrix tests/toy_data/mock_results/flair/counts_matrix.tsv \
    --cantatadb tests/toy_data/annotations/cantata_test.bed \
    --te-bed tests/toy_data/annotations/te_test.bed \
    --min-consensus 2 --max-orf-aa 100 --cpat-threshold 0.39 \
    --out-gtf results_local/lncrna_demo_final.gtf \
    --out-bed results_local/lncrna_demo_final.bed \
    --out-counts results_local/lncrna_demo_counts.tsv \
    --out-report results_local/lncrna_demo_report.tsv
```

### A.2 Ejecutar Epitx Consensus (demo con mock data)

```bash
python scripts/epitx_consensus.py \
    --eligos-dir tests/toy_data/mock_results/eligos \
    --m6anet-dir tests/toy_data/mock_results/m6anet \
    --min-tools 1 --pval 0.05 \
    --out-consensus results_local/epitx_demo_consensus.tsv \
    --out-biotype results_local/epitx_demo_biotype.tsv
```

### A.3 Correr Suite Completa de Tests

```bash
python -m pytest tests/ -v
# Esperado: 136 passed
```

---

## MODO B: Pipeline Completo via Docker en tu PC (con reference data)

Necesitas descargar el reference genome y annotations.

### B.1 Crear estructura de directorios

```bash
cd /c/Users/geope/Downloads/K-CHOPORE-master\ \(2\)/K-CHOPORE-master

mkdir -p data/reference/genome
mkdir -p data/reference/annotations
mkdir -p data/reference/transcriptome
```

### B.2 Descargar TAIR10 reference

```bash
# Genoma TAIR10
cd data/reference/genome
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas
# Renombrar
mv TAIR10_chr_all.fas TAIR10_chr_all.fas.fasta

cd ../../..
```

### B.3 Descargar anotaciones

```bash
# AtRTDv2 GTF (necesario para alignment y FLAIR)
# Descargar manualmente de: https://ics.hutton.ac.uk/atRTD/
# y colocar en data/reference/annotations/AtRTDv2_QUASI_19April2016.gtf

# Anotaciones v3 (Araport11, CANTATAdb, miRBase, PmiREN, TEs)
bash scripts/download_annotations.sh
```

### B.4 Copiar datos ya procesados

Si tienes resultados del servidor, copialos:

```bash
# Copiar FLAIR counts (ya los tienes en results_local/)
mkdir -p results/flair
cp results_local/flair/counts_matrix.tsv results/flair/
cp results_local/flair/reads_manifest.tsv results/flair/

# Copiar DESeq2
mkdir -p results/deseq2
cp results_local/deseq2/* results/deseq2/
```

### B.5 Dry run (ver que haria Snakemake)

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
    -v "$(pwd)":/workspace \
    -w /workspace \
    k-chopore:latest \
    snakemake -n --printshellcmds --configfile config/config.yml 2>&1 | head -50
```

### B.6 Ejecutar pipeline completo

```bash
MSYS_NO_PATHCONV=1 docker run --rm \
    -v "$(pwd)":/workspace \
    -w /workspace \
    k-chopore:latest \
    snakemake --cores 8 --latency-wait 60 --printshellcmds --keep-going \
    --configfile config/config.yml
```

### B.7 Ejecutar solo un modulo

```bash
# Solo lncRNA discovery (Module 1)
MSYS_NO_PATHCONV=1 docker run --rm \
    -v "$(pwd)":/workspace \
    -w /workspace \
    k-chopore:latest \
    snakemake --cores 8 --latency-wait 60 \
    results/lncrna/lncrna_final.gtf

# Solo epitx consensus (Module 4)
MSYS_NO_PATHCONV=1 docker run --rm \
    -v "$(pwd)":/workspace \
    -w /workspace \
    k-chopore:latest \
    snakemake --cores 8 --latency-wait 60 \
    results/epitx/modification_consensus.tsv
```

---

## MODO C: Pipeline en Servidor Remoto (donde estan los BAMs y FAST5)

Este es el modo recomendado para produccion.

### C.1 Conectar al servidor

```bash
ssh usuario@servidor
cd /ruta/a/K-CHOPORE
```

### C.2 Actualizar codigo

```bash
git pull origin main
```

### C.3 Reconstruir imagen Docker (necesario para v3 tools)

```bash
sudo docker build -t k-chopore:v3 .
# Tarda ~30-45 minutos
```

### C.4 Ejecutar Fase A (lncRNA + Epitx) con datos DRS existentes

Verificar que config.yml tiene:
```yaml
params:
  run_lncrna: true
  run_epitx_enhanced: true
  run_smallrna: false
  run_mirna_targets: false
  run_integration: false
```

Ejecutar:
```bash
sudo docker run --rm \
    -v "$(pwd)":/workspace \
    -w /workspace \
    k-chopore:v3 \
    snakemake --cores 40 --latency-wait 60 --printshellcmds --keep-going \
    results/lncrna/lncrna_final.gtf \
    results/lncrna/deseq2/deseq2_genotype_results.csv \
    results/epitx/modification_consensus.tsv
```

### C.5 Ejecutar Fase B (sRNA + Targets) cuando tengas datos Illumina

Primero configura las muestras de sRNA-seq en config.yml:
```yaml
srna_samples:
  srna_wt_c_r1:
    genotype: "WT"
    treatment: "C"
    fastq: "data/raw/srna/srna_wt_c_r1.fastq.gz"
  # ... mas muestras

degradome_samples:
  deg_wt:
    genotype: "WT"
    fastq: "data/raw/degradome/degradome_wt.fastq.gz"

params:
  run_smallrna: true
  run_mirna_targets: true
```

```bash
sudo docker run --rm \
    -v "$(pwd)":/workspace \
    -w /workspace \
    k-chopore:v3 \
    snakemake --cores 40 --latency-wait 60 --printshellcmds --keep-going \
    results/smallrna/shortstack/Results.txt \
    results/targets/target_evidence_table.tsv
```

### C.6 Ejecutar Fase C (Integracion) despues de A+B

```yaml
params:
  run_integration: true
```

```bash
sudo docker run --rm \
    -v "$(pwd)":/workspace \
    -w /workspace \
    k-chopore:v3 \
    snakemake --cores 40 --latency-wait 60 --printshellcmds --keep-going \
    results/integration/wgcna_modules.tsv
```

---

## Problemas Comunes

| Problema | Solucion |
|----------|----------|
| `MSYS_NO_PATHCONV` error | En Git Bash: `export MSYS_NO_PATHCONV=1` antes de docker |
| Docker "permission denied" | Usa `sudo` en Linux, o anade tu usuario al grupo docker |
| Snakemake "MissingInputException" | Faltan archivos de input. Usa `-n` (dry run) para ver que necesita |
| ELIGOS2 "CMH test failure" | Bug conocido de rpy2/R. El pipeline crea placeholders |
| FLAIR "underscore in ID" | Automatico: la pipeline convierte `_` a `-` en el manifest |
| "FileExistsError .snakemake" | Pre-crear: `mkdir -p .snakemake/locks .snakemake/metadata` |
| Docker lento en NTFS | Usa `--latency-wait 60`. Considera WSL2 para mejor rendimiento |
