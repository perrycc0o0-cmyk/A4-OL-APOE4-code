# Software environment

The archived CFFF analysis session records the following tested platform:

- CFFF Linux server;
- Ubuntu 22.04.3 LTS, x86_64;
- R 4.5.2;
- 8 CPU cores and 250 GB RAM for the full mouse preprocessing workflow.

The portable path configuration and OS-aware parallel branches also support
Windows 10/11. The lightweight demo is suitable for a normal desktop. Full
Seurat processing and SCENIC were run on the CFFF Linux server; a complete
end-to-end Windows rerun has not been claimed.

`cfff_linux_session_info.txt` is the available session record from the same
CFFF analysis environment. It includes exact versions for the core Seurat and
tidyverse packages loaded in that run. `package_requirements.tsv` lists every
direct package used in the public scripts and marks whether its exact tested
version was captured. Versions marked `not captured` must be filled from a
fresh server run before the dependency-version checklist item is considered
complete.

Install the declared dependencies with:

```bash
Rscript environment/install_dependencies.R
```

On a typical broadband connection, installation usually takes approximately
15-45 minutes, excluding installation of R itself and download of the large
species-specific cisTarget databases. This is an estimate, not a benchmark.
Run the following on the CFFF server after installation to generate a complete
version record:

```bash
Rscript environment/capture_session_info.R
```

The full workflow requires substantial memory. The archived run used 250 GB;
SCENIC runtime and memory depend strongly on the number of sampled cells and
genes. The code does not require a GPU.
