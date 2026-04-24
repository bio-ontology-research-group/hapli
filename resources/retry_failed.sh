#!/usr/bin/env bash
# Second pass for papers whose primary source 403'd. Tries PMC, arXiv, and publisher mirrors.
set -u
cd "$(dirname "$0")"
LOG="download_papers.log"
UA='Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36'

fetch() {
    local url="$1" out="$2"
    if [ -s "$out" ]; then
        echo "SKIP  $out" >> "$LOG"
        return 0
    fi
    mkdir -p "$(dirname "$out")"
    if curl -fsSL -A "$UA" --max-time 60 -o "$out.tmp" "$url" 2>>"$LOG"; then
        local sz=$(stat -c%s "$out.tmp" 2>/dev/null || echo 0)
        if [ "$sz" -lt 20000 ] && head -c 200 "$out.tmp" | grep -qi '<html\|<!doctype'; then
            echo "HTML  $out ($sz bytes; deleted)" >> "$LOG"
            rm -f "$out.tmp"
            return 1
        fi
        mv "$out.tmp" "$out"
        echo "OK    $out ($sz bytes, retry)" >> "$LOG"
        return 0
    else
        rm -f "$out.tmp"
        echo "FAIL-retry  $out <- $url" >> "$LOG"
        return 1
    fi
}

# PMC mirrors for OUP papers
fetch 'https://europepmc.org/article/MED/28205675/PDF'                                                   'sota/BCFtools-csq_Danecek_2017_Bioinformatics.pdf'
fetch 'https://europepmc.org/article/MED/29750242/PDF'                                                   'algorithms/minimap2_Li_2018_Bioinformatics.pdf'
fetch 'https://europepmc.org/article/MED/33325513/PDF'                                                   'algorithms/Liftoff_Shumate_2021_Bioinformatics.pdf'
fetch 'https://europepmc.org/article/MED/36648328/PDF'                                                   'algorithms/miniprot_Li_2023_Bioinformatics.pdf'
fetch 'https://europepmc.org/article/MED/33590861/PDF'                                                   'algorithms/samtools-bcftools_Danecek_2021_GigaScience.pdf'

# PMC direct PDF route
fetch 'https://pmc.ncbi.nlm.nih.gov/articles/PMC6371407/pdf/nihms-1012268.pdf'                           'algorithms/SpliceAI_Jaganathan_2019_Cell.pdf'
fetch 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8935012/pdf/evac014.pdf'                                 'sota/CompensatoryFrameshift_PMC8935012.pdf'
fetch 'https://pmc.ncbi.nlm.nih.gov/articles/PMC7706150/pdf/ocaa164.pdf'                                 'sota/OpenCRAVAT_Pagel_2020_JAMIA.pdf'

# arXiv mirrors for ML papers
fetch 'https://arxiv.org/pdf/1708.01492'                                                                 'algorithms/minimap2_Li_2018_arxiv.pdf'
fetch 'https://arxiv.org/pdf/2205.13760'                                                                 'algorithms/ESM2_Lin_2022_arxiv.pdf'
fetch 'https://arxiv.org/pdf/2106.15110'                                                                 'algorithms/ESM1v_Meier_2021_arxiv.pdf'
fetch 'https://arxiv.org/pdf/2204.10523'                                                                 'algorithms/ESM-IF1_Hsu_2022_arxiv.pdf'
fetch 'https://arxiv.org/pdf/2312.14222'                                                                 'algorithms/ProteinGym_Notin_2023_arxiv.pdf'

# bioRxiv alternative
fetch 'https://www.biorxiv.org/content/10.1101/2023.09.20.558629v2.full.pdf'                             'algorithms/SaProt_Su_2024_biorxiv.pdf'

# Other mirrors
fetch 'https://zhanggroup.org/TM-align/TM-align.pdf'                                                     'algorithms/TM-align_Zhang_2005_NAR.pdf'

# TOGA preprint on biorxiv
fetch 'https://www.biorxiv.org/content/10.1101/2022.09.08.507143v1.full.pdf'                             'algorithms/TOGA_Kirilenko_2023_biorxiv.pdf'

# AlphaMissense — the biorxiv was rejected; try the Science supplementary route or figshare
fetch 'https://storage.googleapis.com/dm_alphamissense/AlphaMissense_paper.pdf'                          'sota/AlphaMissense_Cheng_2023.pdf'

ok=$(grep -c '^OK.*retry)' "$LOG" || true)
fail=$(grep -c '^FAIL-retry' "$LOG" || true)
echo "=== Retry summary: OK=$ok FAIL=$fail ===" | tee -a "$LOG"
