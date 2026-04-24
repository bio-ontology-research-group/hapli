# Reference Paper Manifest

Auto-downloaded by `download_papers.sh` + `retry_failed.sh`. See `download_papers.log` for per-URL fetch status.

## resources/sota/ — prior art / competitors

| Paper | Status | Note |
|---|---|---|
| Haplosaurus (McLaren 2018, Nat Commun) | ✓ downloaded | |
| VEP (McLaren 2016, Genome Biol) | ✓ downloaded | |
| gnomAD constraint / LOFTEE (Karczewski 2020, Nature) | ✓ downloaded | |
| PRESCOTT (Tekpinar 2025, Genome Biol) | ✓ downloaded | |
| LoFTK (2023, BioData Mining) | ✓ downloaded | |
| BCFtools/csq (Danecek 2017, Bioinformatics) | ✗ manual | OUP blocks curl; try https://academic.oup.com/bioinformatics/article/33/13/2037/3000373 from a browser or via institutional proxy |
| AlphaMissense (Cheng 2023, Science) | ✗ manual | Paywalled; author preprint: https://www.biorxiv.org/content/10.1101/2023.08.29.555570v1 |
| OpenCRAVAT (Pagel 2020, JAMIA) | ✗ manual | Try https://pmc.ncbi.nlm.nih.gov/articles/PMC7706150/ |
| Compensatory frameshift (PMC8935012) | ✗ manual | PMC link returned wrapper page; try landing page directly |
| Protein dynamics epistasis (PNAS 2023) | ✗ manual | PNAS paywall; try https://www.pnas.org/doi/full/10.1073/pnas.2313099121 |

## resources/algorithms/ — tools we integrate

| Paper | Status | Note |
|---|---|---|
| minimap2 (Li 2018) | ✓ arXiv version | `minimap2_Li_2018_arxiv.pdf` (1708.01492) |
| hifiasm (Cheng 2021) | ✓ downloaded | |
| Verkko (Rautiainen 2023) | ✓ downloaded | |
| WhatsHap (Martin 2016) | ✓ downloaded | |
| HPRC (Liao 2023) | ✓ downloaded | |
| Minigraph-Cactus (Hickey 2023) | ✓ downloaded | |
| PanGenie (Ebler 2022) | ✓ downloaded | |
| AlphaFold2 (Jumper 2021) | ✓ downloaded | |
| AlphaFold3 (Abramson 2024) | ✓ downloaded | |
| EVE (Frazer 2021) | ✓ downloaded | |
| Enformer (Avsec 2021) | ✓ downloaded | |
| SpliceTransformer (2024) | ✓ downloaded | |
| MaveDB (Esposito 2019) | ✓ downloaded | |
| SaProt (Su 2024, ICLR) | ✓ downloaded (OpenReview) | |
| ESM2 (Lin 2023) | ✓ arXiv version | `ESM2_Lin_2022_arxiv.pdf` (2205.13760) |
| ESM1v (Meier 2021) | ✓ arXiv version | `ESM1v_Meier_2021_arxiv.pdf` (2106.15110) |
| ESM-IF1 (Hsu 2022) | ✓ arXiv version | `ESM-IF1_Hsu_2022_arxiv.pdf` (2204.10523) |
| SpliceAI (Jaganathan 2019, Cell) | ✗ manual | Cell paywall; try https://pmc.ncbi.nlm.nih.gov/articles/PMC6371407/ from a browser |
| TOGA (Kirilenko 2023, Science) | ✗ manual | Science paywall; try author's GitHub or preprint |
| Liftoff (Shumate 2021, Bioinformatics) | ✗ manual | OUP blocks curl; try https://academic.oup.com/bioinformatics/article/37/12/1639/6035128 |
| miniprot (Li 2023, Bioinformatics) | ✗ manual | OUP blocks curl |
| samtools/bcftools (Danecek 2021, GigaScience) | ✗ manual | OUP blocks curl |
| ProteinGym (Notin 2023, NeurIPS) | ✗ manual | OpenReview PDF; try https://proceedings.neurips.cc/ |
| TM-align (Zhang 2005, NAR) | ✗ manual | OUP blocks curl |

## Retrieval tips

- **Oxford University Press (OUP)** papers (Bioinformatics, NAR, GigaScience) block default `curl` User-Agents. Browser download or an authenticated session (e.g. institutional VPN) generally works.
- **Cell** and **Science** are paywalled; author-hosted preprints or PMC archives are the workaround where available.
- Papers we could not programmatically retrieve are listed in `download_papers.log` with full failed URLs for easy re-attempt.
