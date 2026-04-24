#!/usr/bin/env bash
# Download reference papers for hapli into resources/sota/ and resources/algorithms/.
# Idempotent: skips files that already exist and are non-empty.
# Logs successes and failures to resources/download_papers.log.
# Paywalled papers (Science, Cell, some Nature) will 403 — flagged for manual fetch.

set -u
cd "$(dirname "$0")"
LOG="download_papers.log"
: > "$LOG"

UA='Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0 Safari/537.36'

fetch() {
    local url="$1" out="$2"
    if [ -s "$out" ]; then
        echo "SKIP  $out (already present, $(stat -c%s "$out") bytes)" >> "$LOG"
        return 0
    fi
    mkdir -p "$(dirname "$out")"
    if curl -fsSL -A "$UA" --max-time 60 -o "$out.tmp" "$url" 2>>"$LOG"; then
        # Reject obvious HTML paywalls (< 20 KB usually means landing page, not PDF)
        local sz=$(stat -c%s "$out.tmp" 2>/dev/null || echo 0)
        if [ "$sz" -lt 20000 ] && head -c 200 "$out.tmp" | grep -qi '<html\|<!doctype'; then
            echo "HTML  $out ($sz bytes, likely paywall/redirect; deleted)" >> "$LOG"
            rm -f "$out.tmp"
            return 1
        fi
        mv "$out.tmp" "$out"
        echo "OK    $out ($sz bytes)" >> "$LOG"
        return 0
    else
        rm -f "$out.tmp"
        echo "FAIL  $out <- $url" >> "$LOG"
        return 1
    fi
}

# ------------------------------------------------------------------
# resources/sota/ — direct prior art / competitors
# ------------------------------------------------------------------
fetch 'https://www.nature.com/articles/s41467-018-06542-1.pdf'                                           'sota/Haplosaurus_McLaren_2018_NatCommun.pdf'
fetch 'https://genomebiology.biomedcentral.com/counter/pdf/10.1186/s13059-016-0974-4.pdf'                'sota/VEP_McLaren_2016_GenomeBiol.pdf'
fetch 'https://academic.oup.com/bioinformatics/article-pdf/33/13/2037/25155524/btx100.pdf'               'sota/BCFtools-csq_Danecek_2017_Bioinformatics.pdf'
fetch 'https://genomebiology.biomedcentral.com/counter/pdf/10.1186/s13059-025-03581-y.pdf'               'sota/PRESCOTT_Tekpinar_2025_GenomeBiol.pdf'
fetch 'https://www.nature.com/articles/s41586-020-2308-7.pdf'                                            'sota/gnomAD-constraint_Karczewski_2020_Nature.pdf'
fetch 'https://www.biorxiv.org/content/10.1101/2023.08.29.555570v1.full.pdf'                             'sota/AlphaMissense_Cheng_2023_biorxiv.pdf'
fetch 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8935012/pdf/evac014.pdf'                                 'sota/CompensatoryFrameshift_PMC8935012.pdf'
fetch 'https://www.pnas.org/doi/pdf/10.1073/pnas.2313099121'                                             'sota/ProteinDynamicsEpistasis_2023_PNAS.pdf'
fetch 'https://biodatamining.biomedcentral.com/counter/pdf/10.1186/s13040-023-00321-5.pdf'               'sota/LoFTK_2023_BioDataMining.pdf'
fetch 'https://academic.oup.com/jamia/article-pdf/27/11/1794/34224586/ocaa164.pdf'                       'sota/OpenCRAVAT_Pagel_2020_JAMIA.pdf'

# ------------------------------------------------------------------
# resources/algorithms/ — tools we integrate
# ------------------------------------------------------------------
fetch 'https://academic.oup.com/bioinformatics/article-pdf/34/18/3094/25731859/bty191.pdf'               'algorithms/minimap2_Li_2018_Bioinformatics.pdf'
fetch 'https://academic.oup.com/gigascience/article-pdf/10/2/giab008/36332246/giab008.pdf'               'algorithms/samtools-bcftools_Danecek_2021_GigaScience.pdf'
fetch 'https://www.cell.com/cell/pdf/S0092-8674(18)31629-5.pdf'                                          'algorithms/SpliceAI_Jaganathan_2019_Cell.pdf'
fetch 'https://www.science.org/doi/pdf/10.1126/science.abn3107'                                          'algorithms/TOGA_Kirilenko_2023_Science.pdf'
fetch 'https://academic.oup.com/bioinformatics/article-pdf/37/12/1639/39344866/btaa1016.pdf'             'algorithms/Liftoff_Shumate_2021_Bioinformatics.pdf'
fetch 'https://academic.oup.com/bioinformatics/article-pdf/39/1/btad014/49003112/btad014.pdf'            'algorithms/miniprot_Li_2023_Bioinformatics.pdf'
fetch 'https://www.nature.com/articles/s41592-020-01056-5.pdf'                                           'algorithms/hifiasm_Cheng_2021_NatMethods.pdf'
fetch 'https://www.nature.com/articles/s41587-023-01662-6.pdf'                                           'algorithms/Verkko_Rautiainen_2023_NatBiotechnol.pdf'
fetch 'https://genome.cshlp.org/content/27/5/801.full.pdf'                                               'algorithms/WhatsHap_Martin_2016_GenomeRes.pdf'
fetch 'https://www.biorxiv.org/content/10.1101/2021.07.09.450648v1.full.pdf'                             'algorithms/ESM1v_Meier_2021_NeurIPS_biorxiv.pdf'
fetch 'https://www.science.org/doi/pdf/10.1126/science.ade2574'                                          'algorithms/ESM2_Lin_2023_Science.pdf'
fetch 'https://www.biorxiv.org/content/10.1101/2022.04.10.487779v2.full.pdf'                             'algorithms/ESM-IF1_Hsu_2022_biorxiv.pdf'
fetch 'https://openreview.net/pdf?id=6MRm3G4NiU'                                                         'algorithms/SaProt_Su_2024_ICLR.pdf'
fetch 'https://www.nature.com/articles/s41586-021-04043-8.pdf'                                           'algorithms/EVE_Frazer_2021_Nature.pdf'
fetch 'https://www.nature.com/articles/s41586-021-03819-2.pdf'                                           'algorithms/AlphaFold2_Jumper_2021_Nature.pdf'
fetch 'https://www.nature.com/articles/s41586-024-07487-w.pdf'                                           'algorithms/AlphaFold3_Abramson_2024_Nature.pdf'
fetch 'https://www.nature.com/articles/s41586-023-05896-x.pdf'                                           'algorithms/HPRC_Liao_2023_Nature.pdf'
fetch 'https://www.nature.com/articles/s41587-023-01793-w.pdf'                                           'algorithms/Minigraph-Cactus_Hickey_2023_NatBiotechnol.pdf'
fetch 'https://www.nature.com/articles/s41588-022-01043-w.pdf'                                           'algorithms/PanGenie_Ebler_2022_NatGenet.pdf'
fetch 'https://www.biorxiv.org/content/10.1101/2023.12.07.570727v1.full.pdf'                             'algorithms/ProteinGym_Notin_2023_NeurIPS_biorxiv.pdf'
fetch 'https://genomebiology.biomedcentral.com/counter/pdf/10.1186/s13059-019-1845-6.pdf'                'algorithms/MaveDB_Esposito_2019_GenomeBiol.pdf'
fetch 'https://academic.oup.com/nar/article-pdf/33/7/2302/7155478/gki524.pdf'                            'algorithms/TM-align_Zhang_2005_NAR.pdf'
fetch 'https://www.nature.com/articles/s41467-024-53088-6.pdf'                                           'algorithms/SpliceTransformer_2024_NatCommun.pdf'
fetch 'https://www.nature.com/articles/s41592-021-01252-x.pdf'                                           'algorithms/Enformer_Avsec_2021_NatMethods.pdf'

# ------------------------------------------------------------------
# Summary
# ------------------------------------------------------------------
ok=$(grep -c '^OK' "$LOG" || true)
skip=$(grep -c '^SKIP' "$LOG" || true)
fail=$(grep -c '^FAIL\|^HTML' "$LOG" || true)
{
    echo
    echo "=== Summary ==="
    echo "OK:   $ok"
    echo "SKIP: $skip"
    echo "FAIL: $fail"
    echo
    echo "Failures need manual retrieval:"
    grep -E '^FAIL|^HTML' "$LOG" || true
} | tee -a "$LOG"
