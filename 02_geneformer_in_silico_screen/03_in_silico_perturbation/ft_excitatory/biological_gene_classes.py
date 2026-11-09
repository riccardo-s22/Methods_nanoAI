"""
Curated gene-class and artifact-set definitions for the excitatory-neuron
biological-candidate filter (scripts/07_biological_candidate_filter.py).

Class membership is exclusive and resolved in the order classes appear in
GENE_CLASS_SETS. If a symbol matches none, gene_class = "other".

These lists are deliberately conservative — they cover well-known members
of each class likely to appear in a perturbation screen of human cortical
neurons. They are not exhaustive references; edit freely.
"""
from __future__ import annotations

# ── Artifact / contamination filters ──────────────────────────────────────
ARTIFACT_PREFIXES = ("MT-", "RPL", "RPS")          # mitochondrial-encoded, cytoplasmic ribosomal proteins
ARTIFACT_LITERALS = {"MALAT1"}                      # abundant nuclear lncRNA

HEMOGLOBINS = {"HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"}

# Non-neuronal contamination markers — flagged in excitatory-neuron context.
NONNEURONAL_MARKERS = {
    # Astrocyte
    "GFAP", "AQP4", "ALDH1L1", "SLC1A2", "SLC1A3",
    # Microglia / immune
    "PTPRC", "C1QA", "C1QB", "C1QC", "TYROBP", "AIF1", "LST1",
    # Oligodendrocyte / myelin
    "PLP1", "MBP", "MOG", "MAG", "MOBP",
    # Endothelial / pericyte
    "CLDN5", "PECAM1", "VWF", "PDGFRB", "RGS5",
}

# Non-microglia contamination markers — flagged in microglia/macrophage
# screens. Excludes myeloid/microglia markers (PTPRC, C1QA/B/C, TYROBP, AIF1,
# LST1, CSF1R, CX3CR1, TMEM119, P2RY12, ITGAM, CD68, CD163, CD14, MRC1,
# FCGR1A, CCR2, ITGB2) which are the biology of interest.
NON_MICROGLIA_MARKERS = {
    # Excitatory-neuron markers
    "SLC17A7", "SLC17A6", "CAMK2A", "NEUROD2", "NEUROD6",
    "SATB2", "TBR1", "BCL11B", "FEZF2", "CUX1", "CUX2",
    # Inhibitory-neuron markers
    "GAD1", "GAD2", "SLC32A1",
    "PVALB", "SST", "VIP", "NPY", "CALB1", "CALB2", "RELN",
    "DLX1", "DLX2", "DLX5", "DLX6",
    # Pan-neuronal markers
    "RBFOX3", "MAP2", "TUBB3", "ENO2", "STMN2",
    "NEFL", "NEFM", "NEFH",
    "SYN1", "SYP", "SYT1", "SNAP25",
    # Astrocyte markers
    "GFAP", "AQP4", "ALDH1L1", "SLC1A2", "SLC1A3", "S100B", "SOX9",
    "FGFR3", "AGT", "ALDOC", "GJA1", "GLUL", "SLC4A4",
    # Oligodendrocyte lineage (mature + OPC)
    "PLP1", "MBP", "MOG", "MAG", "MOBP",
    "OLIG1", "OLIG2", "SOX10",
    "PDGFRA", "CSPG4",
    # Endothelial / pericyte
    "CLDN5", "PECAM1", "VWF", "PDGFRB", "RGS5",
}

# Non-astrocyte contamination markers — flagged in astrocyte screens.
# Excludes astrocyte markers (GFAP, AQP4, ALDH1L1, SLC1A2, SLC1A3, S100B,
# SOX9, FGFR3, AGT, ALDOC, GJA1, GLUL) which are the biology of interest.
NON_ASTROCYTE_MARKERS = {
    # Excitatory-neuron markers
    "SLC17A7", "SLC17A6", "CAMK2A", "NEUROD2", "NEUROD6",
    "SATB2", "TBR1", "BCL11B", "FEZF2", "CUX1", "CUX2",
    # Inhibitory-neuron markers
    "GAD1", "GAD2", "SLC32A1",
    "PVALB", "SST", "VIP", "NPY", "CALB1", "CALB2", "RELN",
    "DLX1", "DLX2", "DLX5", "DLX6",
    # Pan-neuronal markers
    "RBFOX3", "MAP2", "TUBB3", "ENO2", "STMN2",
    "NEFL", "NEFM", "NEFH",
    "SYN1", "SYP", "SYT1", "SNAP25",
    # Microglia / immune
    "PTPRC", "C1QA", "C1QB", "C1QC", "TYROBP", "AIF1", "LST1",
    "CSF1R", "CX3CR1", "TMEM119", "P2RY12", "ITGAM",
    # Oligodendrocyte lineage (mature + OPC)
    "PLP1", "MBP", "MOG", "MAG", "MOBP",
    "OLIG1", "OLIG2", "SOX10",
    "PDGFRA", "CSPG4",
    # Endothelial / pericyte
    "CLDN5", "PECAM1", "VWF", "PDGFRB", "RGS5",
}

# Non-oligodendrocyte contamination markers — flagged in oligodendrocyte
# screens. Excludes mature-oligo markers (PLP1, MBP, MOG, MAG, MOBP) and
# pan-oligo-lineage TFs (SOX10, OLIG1, OLIG2 - shared with OPCs).
NON_OLIGODENDROCYTE_MARKERS = {
    # Excitatory-neuron markers
    "SLC17A7", "SLC17A6", "CAMK2A", "NEUROD2", "NEUROD6",
    "SATB2", "TBR1", "BCL11B", "FEZF2", "CUX1", "CUX2",
    # Inhibitory-neuron markers
    "GAD1", "GAD2", "SLC32A1",
    "PVALB", "SST", "VIP", "NPY", "CALB1", "CALB2", "RELN",
    "DLX1", "DLX2", "DLX5", "DLX6",
    # Pan-neuronal markers
    "RBFOX3", "MAP2", "TUBB3", "ENO2", "STMN2",
    "NEFL", "NEFM", "NEFH",
    "SYN1", "SYP", "SYT1", "SNAP25",
    # Astrocyte markers
    "GFAP", "AQP4", "ALDH1L1", "SLC1A2", "SLC1A3", "S100B", "SOX9",
    # Microglia / immune
    "PTPRC", "C1QA", "C1QB", "C1QC", "TYROBP", "AIF1", "LST1",
    "CSF1R", "CX3CR1", "TMEM119", "P2RY12", "ITGAM",
    # OPC-distinguishing markers (close lineage but distinct from mature oligo)
    "PDGFRA", "CSPG4",
    # Endothelial / pericyte
    "CLDN5", "PECAM1", "VWF", "PDGFRB", "RGS5",
}

# ── Biological gene classes ───────────────────────────────────────────────
# (resolved in this dict's insertion order; first hit wins)

IMMEDIATE_EARLY_STRESS = {
    "FOS", "FOSB", "FOSL1", "FOSL2",
    "JUN", "JUNB", "JUND",
    "ATF3", "ATF4",
    "EGR1", "EGR2", "EGR3", "EGR4",
    "ARC", "NR4A1", "NR4A2", "NR4A3",
    "IER2", "IER3", "IER5", "BTG2",
    "DUSP1", "DUSP2", "DUSP5", "DUSP6",
    "GADD45B", "GADD45G",  # also apoptosis_stress; classified here first
    "MAFF", "MAFG", "MAFK",
    "CYR61", "CCN1", "CTGF", "CCN2",
    "NPAS4",
}

MECHANOTRANSDUCTION = {
    # Mechanosensitive ion channels (canonical mechanosensors)
    "PIEZO1", "PIEZO2",
    "TRPV2", "TRPV4", "TRPA1",
    "TRPM3", "TRPM7",
    "TRPC1", "TRPC3", "TRPC6",
    "KCNK2", "KCNK4", "KCNK10",          # TREK-1 / TRAAK / TREK-2 — stretch/pressure
    "ASIC1", "ASIC2", "ASIC3", "ASIC4", "ASIC5",
    "TMC1", "TMC2", "TMC4", "TMC6", "TMC7", "TMC8",
    "TMEM150C", "ELKIN1",
    # Hippo / YAP-TAZ mechanosensing
    "YAP1", "WWTR1",                      # WWTR1 = TAZ
    "LATS1", "LATS2",
    "STK3", "STK4",                       # MST2, MST1
    "SAV1", "MOB1A", "MOB1B",
    "NF2",                                # Merlin
    "TEAD1", "TEAD2", "TEAD3", "TEAD4",
    "VGLL1", "VGLL2", "VGLL3", "VGLL4",
    "AMOT", "AMOTL1", "AMOTL2",
    "FAT1", "FAT2", "FAT3", "FAT4",
    "DCHS1", "DCHS2",
    "CRB1", "CRB2", "CRB3",
    "MPP5", "PATJ", "PALS1",
    # SRF / MRTF (mechanically-induced transcription)
    "SRF",
    "MKL1", "MKL2", "MRTFA", "MRTFB",
    # Integrins
    "ITGA1", "ITGA2", "ITGA3", "ITGA4", "ITGA5", "ITGA6", "ITGA7", "ITGA8",
    "ITGA9", "ITGA10", "ITGA11", "ITGAV", "ITGAM", "ITGAL", "ITGAX", "ITGAE",
    "ITGB1", "ITGB2", "ITGB3", "ITGB4", "ITGB5", "ITGB6", "ITGB7", "ITGB8",
    # Focal adhesion / mechanosensitive adaptors
    "PTK2", "PTK2B",                      # FAK / PYK2
    "PXN", "LPP",
    "VCL", "ZYX",
    "TNS1", "TNS2", "TNS3", "TNS4",
    "TLN1", "TLN2",
    "ILK", "LIMS1", "LIMS2", "PARVA", "PARVB", "PARVG", "RSU1",
    "CRK", "CRKL", "DOCK1",
    "BCAR1", "BCAR3",
    # Dystroglycan complex — couples ECM (laminin) to cytoskeleton (mechano-bridge)
    "DAG1", "DMD", "UTRN",
    "SGCA", "SGCB", "SGCG", "SGCD",
    # Actin / actin-binding (cytoskeletal mechanics)
    "ACTA1", "ACTA2", "ACTB", "ACTC1", "ACTG1", "ACTG2",
    "ACTN1", "ACTN2", "ACTN3", "ACTN4",
    "FLNA", "FLNB", "FLNC",
    "CFL1", "CFL2", "DSTN",               # ADF/cofilin family
    "GSN", "VILL",
    "EZR", "RDX", "MSN",                  # ERM scaffolds
    "TLN2",
    # Myosin / contractile machinery
    "MYL9", "MYL12A", "MYL12B",
    "MYH9", "MYH10", "MYH11", "MYH14",
    "MYLK", "MYLK2", "MYLK3", "MYLK4",
    # Formins, nucleation
    "DIAPH1", "DIAPH2", "DIAPH3",
    "FMN1", "FMN2",
    "ARPC1A", "ARPC1B", "ARPC2", "ARPC4", "ARPC5", "ARPC5L",  # ARPC3 omitted (overlap with IEG-like)
    "WASF1", "WASF2", "WASF3",
    "WASL", "WAS",
    "ACTR2", "ACTR3", "ACTR3B", "ACTR3C",
    # Rho-GTPase signaling
    "RHOA", "RHOB", "RHOC",
    "ROCK1", "ROCK2",
    "RAC1", "RAC2", "RAC3",
    "CDC42", "CDC42BPA", "CDC42BPB", "CDC42BPG",
    "ARHGAP1", "ARHGAP4", "ARHGAP6", "ARHGAP17", "ARHGAP18", "ARHGAP21",
    "ARHGAP24", "ARHGAP26", "ARHGAP29", "ARHGAP32", "ARHGAP35", "ARHGAP44",
    "ARHGEF1", "ARHGEF2", "ARHGEF7", "ARHGEF10", "ARHGEF11", "ARHGEF12",
    "ARHGEF17", "ARHGEF18", "ARHGEF28",
    "LARG", "VAV1", "VAV2", "VAV3",
    # Membrane repair / plasma-membrane mechanics
    "AHNAK", "AHNAK2",                    # mechanically-induced membrane scaffolds
    "DYSF",                                # dysferlin — sarcolemmal repair
    "TRIM72", "MG53",
    "ANXA1", "ANXA2", "ANXA5", "ANXA6",
    "CAV1", "CAV2", "CAV3",                # caveolae — mechanoprotection
    "CAVIN1", "CAVIN2", "CAVIN3", "CAVIN4",
    "PTRF",
    # CCN family (Hippo/YAP downstream, mechanical-stress induced)
    # NB: CCN1/CYR61 and CCN2/CTGF stay in IEG (earlier in resolution order);
    # CCN3-6 land here.
    "CCN3", "NOV",
    "CCN4", "WISP1",
    "CCN5", "WISP2",
    "CCN6", "WISP3",
    # Nuclear mechanotransduction (lamins, LINC complex)
    "LMNA", "LMNB1", "LMNB2",
    "SYNE1", "SYNE2", "SYNE3", "SYNE4",
    "SUN1", "SUN2",
    "EMD",
    "TMEM43",
    # Mechanical-stress kinases / signaling
    "MAP3K8", "MAP4K4",
    "JNK1", "JNK2", "JNK3", "MAPK8", "MAPK9", "MAPK10",   # JNK family — stress + stretch
    "MAPK11", "MAPK12", "MAPK13", "MAPK14",                # p38 — stretch-activated
}

ECM = {
    # ── Collagens (all 28 types) ────────────────────────────────────────
    "COL1A1", "COL1A2",
    "COL2A1",
    "COL3A1",
    "COL4A1", "COL4A2", "COL4A3", "COL4A4", "COL4A5", "COL4A6",
    "COL5A1", "COL5A2", "COL5A3",
    "COL6A1", "COL6A2", "COL6A3", "COL6A5", "COL6A6",
    "COL7A1",
    "COL8A1", "COL8A2",
    "COL9A1", "COL9A2", "COL9A3",
    "COL10A1",
    "COL11A1", "COL11A2",
    "COL12A1",
    "COL13A1", "COL14A1", "COL15A1", "COL16A1", "COL17A1", "COL18A1",
    "COL19A1", "COL20A1", "COL21A1", "COL22A1", "COL23A1", "COL24A1",
    "COL25A1", "COL26A1", "COL27A1", "COL28A1",
    # ── Laminins ───────────────────────────────────────────────────────
    "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5",
    "LAMB1", "LAMB2", "LAMB3", "LAMB4",
    "LAMC1", "LAMC2", "LAMC3",
    # ── Basement-membrane scaffolds ─────────────────────────────────────
    "NID1", "NID2",
    "HSPG2",                      # perlecan
    "AGRN",                       # agrin proteoglycan
    # ── Fibronectin & fibrillin family ──────────────────────────────────
    "FN1",
    "FBN1", "FBN2", "FBN3",
    "FBLN1", "FBLN2", "EFEMP1", "FBLN5", "EFEMP2", "FBLN7",
    "HMCN1", "HMCN2",
    "MFAP1", "MFAP2", "MFAP3", "MFAP4", "MFAP5",
    # ── Elastin ────────────────────────────────────────────────────────
    "ELN",
    "EMILIN1", "EMILIN2", "EMILIN3",
    "MFGE8",
    # ── Heparan-sulfate proteoglycans (cell-surface) ────────────────────
    "SDC1", "SDC2", "SDC3", "SDC4",
    "GPC1", "GPC2", "GPC3", "GPC4", "GPC5", "GPC6",
    # ── Chondroitin-sulfate proteoglycans (PNN core) ────────────────────
    "ACAN", "VCAN", "BCAN", "NCAN", "CSPG4", "CSPG5",
    "PHLDA1",
    # ── Hyaluronan: synthases, degradation, link proteins, receptors ────
    "HAS1", "HAS2", "HAS3",
    "HYAL1", "HYAL2", "HYAL3", "HYAL4",
    "HAPLN1", "HAPLN2", "HAPLN3", "HAPLN4",
    "CD44",
    "STAB1", "STAB2", "LYVE1",
    # ── Small leucine-rich proteoglycans (SLRPs) ────────────────────────
    "DCN", "BGN", "FMOD", "LUM", "ASPN", "ECM2",
    "OGN", "OMD", "PRELP", "KERA", "OPTC", "CHAD", "NYX",
    # ── Matricellular / matricrine proteins ────────────────────────────
    "SPARC", "SPARCL1",                          # SPARCL1 = hevin
    "THBS1", "THBS2", "THBS3", "THBS4",          # thrombospondins
    "TNC", "TNN", "TNR", "TNXB",                 # tenascins (TNR = brain-enriched)
    "POSTN",
    "MGP", "MEPE", "DMP1", "SPP1",               # SIBLINGs / osteopontin
    "CTHRC1", "COMP",
    "OLFM1", "OLFM2", "OLFM3", "OLFM4",          # olfactomedins (noelins)
    "OLFML1", "OLFML2A", "OLFML2B", "OLFML3",
    # ── MMPs, TIMPs, ADAM/ADAMTS proteases ──────────────────────────────
    "MMP1", "MMP2", "MMP3", "MMP7", "MMP8", "MMP9", "MMP10",
    "MMP11", "MMP12", "MMP13", "MMP14", "MMP15", "MMP16", "MMP17",
    "MMP19", "MMP20", "MMP21", "MMP23A", "MMP23B",
    "MMP24", "MMP25", "MMP26", "MMP27", "MMP28",
    "TIMP1", "TIMP2", "TIMP3", "TIMP4",
    "ADAMTS1", "ADAMTS2", "ADAMTS3", "ADAMTS4", "ADAMTS5", "ADAMTS6",
    "ADAMTS7", "ADAMTS8", "ADAMTS9", "ADAMTS10", "ADAMTS12", "ADAMTS13",
    "ADAMTS14", "ADAMTS15", "ADAMTS16", "ADAMTS17", "ADAMTS18",
    "ADAMTS19", "ADAMTS20",
    "ADAMTSL1", "ADAMTSL2", "ADAMTSL3", "ADAMTSL4", "ADAMTSL5",
    "ADAM8", "ADAM9", "ADAM10", "ADAM11", "ADAM12", "ADAM15", "ADAM17",
    "ADAM19", "ADAM21", "ADAM22", "ADAM23", "ADAM28", "ADAM29", "ADAM33",
    # ── ECM cross-linking & post-translational modifiers ────────────────
    "LOX", "LOXL1", "LOXL2", "LOXL3", "LOXL4",
    "TGM1", "TGM2", "TGM3", "TGM4", "TGM5", "TGM6", "TGM7",
    "PLOD1", "PLOD2", "PLOD3",
    "P4HA1", "P4HA2", "P4HA3", "P4HB", "P4HTM",
    "PLG", "PLAT", "PLAU", "PLAUR", "SERPINE1",   # plasmin axis (ECM turnover)
    "SERPINH1",                                    # HSP47 — collagen chaperone
    # ── Brain ECM / Reelin signaling ───────────────────────────────────
    "RELN", "DAB1", "VLDLR", "LRP8",
    "MATN1", "MATN2", "MATN3", "MATN4",            # matrilins
    # ── Other ECM-related ──────────────────────────────────────────────
    "MMRN1", "MMRN2",
    "EMID1", "EMID2",
    "VWF",   # endothelial, also in NONNEURONAL_MARKERS — artifact filter removes first
}

SYNAPTIC_PLASTICITY = {
    "BDNF", "NTRK2", "NGF", "NTRK1",
    "HOMER1", "HOMER2", "HOMER3",
    "DLG1", "DLG2", "DLG3", "DLG4",   # PSD95 family
    "SYN1", "SYN2", "SYN3",
    "SYT1", "SYT2", "SYT4", "SYT5", "SYT7", "SYT9", "SYT11",
    "SYP", "SYPL1", "SYPL2",
    "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B",
    "GRIA1", "GRIA2", "GRIA3", "GRIA4",
    "GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5",
    "GRM1", "GRM2", "GRM3", "GRM4", "GRM5", "GRM7", "GRM8",
    "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",
    "SHANK1", "SHANK2", "SHANK3",
    "NLGN1", "NLGN2", "NLGN3", "NLGN4X", "NLGN4Y",
    "NRXN1", "NRXN2", "NRXN3",
    "PSD95", "PCLO", "BSN", "RIMS1", "RIMS2",
    "SNAP25", "STX1A", "STX1B", "STXBP1", "VAMP1", "VAMP2", "VAMP3",
    "DLGAP1", "DLGAP2", "DLGAP3", "DLGAP4",
    "STX2", "RAB3A", "RAB3B", "RAB3C", "MUNC13", "UNC13A", "UNC13B",
    "DLGAP", "GEPHYRIN", "GPHN",
    "SYNGAP1", "DCC", "ROBO1", "ROBO2",
}

CALCIUM_ION_CHANNEL = {
    # Voltage-gated calcium (alpha, beta, gamma, alpha2/delta)
    "CACNA1A", "CACNA1B", "CACNA1C", "CACNA1D", "CACNA1E", "CACNA1F",
    "CACNA1G", "CACNA1H", "CACNA1I", "CACNA1S",
    "CACNA2D1", "CACNA2D2", "CACNA2D3", "CACNA2D4",
    "CACNB1", "CACNB2", "CACNB3", "CACNB4",
    "CACNG1", "CACNG2", "CACNG3", "CACNG4", "CACNG5", "CACNG6", "CACNG7", "CACNG8",
    # Voltage-gated potassium
    "KCNA1", "KCNA2", "KCNA3", "KCNA4", "KCNA5", "KCNA6", "KCNA7", "KCNA10",
    "KCNB1", "KCNB2", "KCNC1", "KCNC2", "KCNC3", "KCNC4",
    "KCND1", "KCND2", "KCND3",
    "KCNQ1", "KCNQ2", "KCNQ3", "KCNQ4", "KCNQ5",
    "KCNH1", "KCNH2", "KCNH3", "KCNH4", "KCNH5", "KCNH6", "KCNH7", "KCNH8",
    "KCNJ1", "KCNJ2", "KCNJ3", "KCNJ4", "KCNJ5", "KCNJ6", "KCNJ8", "KCNJ10", "KCNJ11", "KCNJ12", "KCNJ15",
    "KCNK1", "KCNK2", "KCNK3", "KCNK4", "KCNK5", "KCNK9", "KCNK10", "KCNK13", "KCNK15",
    "KCNN1", "KCNN2", "KCNN3", "KCNN4",
    "KCNT1", "KCNT2",
    # Voltage-gated sodium
    "SCN1A", "SCN2A", "SCN3A", "SCN4A", "SCN5A", "SCN8A", "SCN9A", "SCN10A", "SCN11A",
    "SCN1B", "SCN2B", "SCN3B", "SCN4B",
    # HCN
    "HCN1", "HCN2", "HCN3", "HCN4",
    # Calmodulin, IP3R, RyR, TRP, calcium-binding
    "CALM1", "CALM2", "CALM3",
    "ITPR1", "ITPR2", "ITPR3",
    "RYR1", "RYR2", "RYR3",
    "TRPC1", "TRPC3", "TRPC4", "TRPC5", "TRPC6", "TRPC7",
    "TRPM2", "TRPM4", "TRPM7",
    "TRPV1", "TRPV2", "TRPV4",
    "STIM1", "STIM2", "ORAI1", "ORAI2", "ORAI3",
    "PMCA", "ATP2B1", "ATP2B2", "ATP2B3", "ATP2B4",
    "SLC8A1", "SLC8A2", "SLC8A3",
    # Voltage-gated chloride
    "CLCN1", "CLCN2", "CLCN3", "CLCN4", "CLCN5", "CLCN6", "CLCN7",
}

TDP43_RNA_SPLICING = {
    "TARDBP", "FUS", "TAF15", "EWSR1",
    "HNRNPA1", "HNRNPA2B1", "HNRNPA3", "HNRNPC", "HNRNPD", "HNRNPF",
    "HNRNPH1", "HNRNPH2", "HNRNPH3", "HNRNPK", "HNRNPL", "HNRNPM", "HNRNPR", "HNRNPU", "HNRNPUL1", "HNRNPUL2",
    "SFPQ", "NONO", "PSPC1",
    "MATR3", "TIA1", "TIAL1",
    "SRSF1", "SRSF2", "SRSF3", "SRSF4", "SRSF5", "SRSF6", "SRSF7", "SRSF9", "SRSF10", "SRSF11", "SRSF12",
    "U2AF1", "U2AF2", "U2AF1L4",
    "SNRPA", "SNRPA1", "SNRPB", "SNRPB2", "SNRPC", "SNRPD1", "SNRPD2", "SNRPD3", "SNRPE", "SNRPF", "SNRPG", "SNRPN",
    "SNRNP70", "SNRNP200", "SNRNP25", "SNRNP35", "SNRNP40", "SNRNP48",
    "PTBP1", "PTBP2", "PTBP3",
    "RBFOX1", "RBFOX2", "RBFOX3",
    "NOVA1", "NOVA2",
    "CELF1", "CELF2", "CELF3", "CELF4", "CELF5", "CELF6",
    "ELAVL1", "ELAVL2", "ELAVL3", "ELAVL4",
    "STAU1", "STAU2",
    "DHX9", "DDX5", "DDX17",
    "TRA2A", "TRA2B",
    "SRRM1", "SRRM2", "SRRM3", "SRRM4",
    "EIF4A3", "MAGOH", "RBM8A",  # EJC
    "UPF1", "UPF2", "UPF3A", "UPF3B",  # NMD-adjacent
    "PRPF8", "PRPF31", "PRPF38A", "PRPF38B", "PRPF39", "PRPF40A", "PRPF40B", "PRPF4B",
    "PCBP1", "PCBP2", "PCBP3", "PCBP4",
    "ATXN2", "ATXN2L",
}

# Nuclear-encoded mitochondrial genes — explicit set + prefix scan.
MITO_NUCLEAR_PREFIXES = (
    "NDUF",          # Complex I subunits / assembly
    "SDHA", "SDHB", "SDHC", "SDHD", "SDHAF",  # Complex II
    "UQCR",          # Complex III
    "COX",           # Complex IV (excluded: COX would also match cyclooxygenase PTGS aliases; PTGS gene symbols are PTGS, not COX)
    "ATP5", "ATPAF", "ATPSCKMT",  # Complex V
    "MRPS", "MRPL",  # mitochondrial ribosomes
    "TOMM", "TIMM",  # import machinery
    "TFAM", "TFB1M", "TFB2M",
    "POLG", "POLG2", "POLRMT",
    "SLC25A",        # mitochondrial carriers
    "BCS1L", "SURF1", "LRPPRC", "SCO1", "SCO2",
    "PINK1", "PARK7", "PRKN",
    "OPA1", "MFN1", "MFN2", "DRP1", "DNM1L",
    "PPARGC1A", "PPARGC1B",
    "NRF1", "NRF2", "GABPA", "GABPB1",
    "FIS1", "MFF",
    "COA",           # cytochrome c oxidase assembly factors COA1..7
)

MITO_NUCLEAR_LITERALS = {
    "CYC1", "CYCS", "ETFA", "ETFB", "ETFDH",
    "HADHA", "HADHB", "ECHS1", "ACADM", "ACADL", "ACADVL", "ACADS",
    "MT-ATP6", "MT-ATP8",  # these are encoded on mtDNA but start with MT- and are filtered;
                              # included here so reclassification wouldn't break if MT- filter removed.
}

PROTEOSTASIS = {
    # HSPs and cochaperones
    "HSPA1A", "HSPA1B", "HSPA2", "HSPA4", "HSPA4L", "HSPA5", "HSPA6", "HSPA7", "HSPA8", "HSPA9", "HSPA12A", "HSPA12B", "HSPA13", "HSPA14",
    "HSPB1", "HSPB2", "HSPB3", "HSPB6", "HSPB7", "HSPB8", "HSPB9", "HSPB11",
    "HSP90AA1", "HSP90AB1", "HSP90B1",
    "HSPE1", "HSPD1",
    "DNAJA1", "DNAJA2", "DNAJA3", "DNAJA4",
    "DNAJB1", "DNAJB2", "DNAJB4", "DNAJB6", "DNAJB9", "DNAJB11", "DNAJB12", "DNAJB14",
    "DNAJC1", "DNAJC2", "DNAJC3", "DNAJC5", "DNAJC6", "DNAJC7", "DNAJC8", "DNAJC10", "DNAJC11", "DNAJC12", "DNAJC13",
    "BAG1", "BAG2", "BAG3", "BAG4", "BAG5", "BAG6",
    # UPS - ubiquitin / ligases / proteasome
    "UBA1", "UBA52", "UBB", "UBC", "UBD",
    "UBE2A", "UBE2B", "UBE2D1", "UBE2D2", "UBE2D3", "UBE2E1", "UBE2E2", "UBE2E3", "UBE2G1", "UBE2G2",
    "UBE2H", "UBE2I", "UBE2J1", "UBE2J2", "UBE2K", "UBE2L3", "UBE2L6", "UBE2M", "UBE2N", "UBE2O", "UBE2Q1", "UBE2R2", "UBE2S", "UBE2T", "UBE2V1", "UBE2V2", "UBE2Z",
    "VCP", "STUB1",
    "PSMA1", "PSMA2", "PSMA3", "PSMA4", "PSMA5", "PSMA6", "PSMA7", "PSMA8",
    "PSMB1", "PSMB2", "PSMB3", "PSMB4", "PSMB5", "PSMB6", "PSMB7", "PSMB8", "PSMB9", "PSMB10", "PSMB11",
    "PSMC1", "PSMC2", "PSMC3", "PSMC4", "PSMC5", "PSMC6",
    "PSMD1", "PSMD2", "PSMD3", "PSMD4", "PSMD5", "PSMD6", "PSMD7", "PSMD8", "PSMD9", "PSMD10", "PSMD11", "PSMD12", "PSMD13", "PSMD14",
    "PSME1", "PSME2", "PSME3", "PSME4",
    "PSMG1", "PSMG2", "PSMG3", "PSMG4",
    # Autophagy
    "ATG3", "ATG4A", "ATG4B", "ATG4C", "ATG4D", "ATG5", "ATG7", "ATG9A", "ATG9B", "ATG10", "ATG12", "ATG13", "ATG14", "ATG16L1", "ATG16L2",
    "MAP1LC3A", "MAP1LC3B", "MAP1LC3B2", "MAP1LC3C", "GABARAP", "GABARAPL1", "GABARAPL2",
    "SQSTM1", "BECN1", "BECN2", "NBR1", "OPTN", "TBK1", "NDP52", "CALCOCO2", "TAX1BP1",
    # UPR
    "ATF6", "ATF6B", "XBP1", "DDIT3", "EIF2AK1", "EIF2AK2", "EIF2AK3", "EIF2AK4", "ERN1", "ERN2",
    "EDEM1", "EDEM2", "EDEM3", "HRD1", "SYVN1", "HERPUD1",
    # Other PQC
    "CCT2", "CCT3", "CCT4", "CCT5", "CCT6A", "CCT6B", "CCT7", "CCT8", "TCP1",
    "PFDN1", "PFDN2", "PFDN3", "PFDN4", "PFDN5", "PFDN6",
}

APOPTOSIS_STRESS = {
    "BAX", "BAK1", "BAD", "BID", "BIK", "BBC3", "PMAIP1", "BMF",
    "BCL2", "BCL2L1", "BCL2L2", "BCL2L11", "MCL1", "BCL2A1",
    "CASP1", "CASP2", "CASP3", "CASP4", "CASP5", "CASP6", "CASP7", "CASP8", "CASP9", "CASP10", "CASP14",
    "TP53", "MDM2", "MDM4", "TP53I3", "TP53BP1", "TP53BP2", "TP63", "TP73",
    "CDKN1A", "CDKN1B", "CDKN1C", "CDKN2A", "CDKN2B", "CDKN2C", "CDKN2D",
    "GADD45A",
    "PARP1", "PARP2", "PARP3", "PARP4",
    "AIFM1", "AIFM2", "AIFM3", "APAF1",
    "XIAP", "BIRC2", "BIRC3", "BIRC5", "BIRC6", "BIRC7", "BIRC8",
    "FAS", "FASLG", "FADD", "TRADD", "TRAF1", "TRAF2", "TRAF3", "TRAF6",
    "RIPK1", "RIPK3", "MLKL",
    "DIABLO", "HTRA2", "ENDOG",
    "AKT1", "AKT2", "AKT3", "PTEN",
    "FOXO1", "FOXO3", "FOXO4", "FOXO6",
    "TNF", "TNFRSF1A", "TNFRSF1B", "TNFRSF10A", "TNFRSF10B", "TNFRSF10C", "TNFRSF10D",
    "BCAP31",
}

TRANSCRIPTION_FACTOR = {
    # Cortical-neuronal TFs
    "TBR1", "TBR2", "EOMES",
    "NEUROD1", "NEUROD2", "NEUROD4", "NEUROD6",
    "POU3F2", "POU3F3", "POU3F4",
    "SATB1", "SATB2",
    "BCL11A", "BCL11B", "CUX1", "CUX2",
    "FOXP1", "FOXP2", "FOXP4", "FOXG1",
    "MEF2A", "MEF2B", "MEF2C", "MEF2D",
    "CREB1", "CREB3", "CREB3L1", "CREB3L2", "CREB3L3", "CREB3L4", "CREBL2", "CREBZF", "CREBBP",
    "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6",
    "NFKB1", "NFKB2", "RELA", "RELB", "REL", "NFKBIA", "NFKBIB", "NFKBIE",
    "E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8",
    "HIF1A", "ARNT", "EPAS1", "HIF3A",
    "SP1", "SP2", "SP3", "SP4", "SP5", "SP6", "SP7", "SP8", "SP9",
    "KLF1", "KLF2", "KLF3", "KLF4", "KLF5", "KLF6", "KLF7", "KLF8", "KLF9", "KLF10", "KLF11", "KLF12", "KLF13", "KLF14", "KLF15", "KLF16", "KLF17",
    "NFIA", "NFIB", "NFIC", "NFIX",
    "NFE2", "NFE2L1", "NFE2L2", "NFE2L3",
    "MYC", "MYCN", "MYCL", "MAX", "MXD1", "MXD3", "MXD4", "MXI1", "MAZ", "MNT", "MGA",
    "DLX1", "DLX2", "DLX3", "DLX4", "DLX5", "DLX6",
    "LHX1", "LHX2", "LHX3", "LHX5", "LHX6", "LHX8", "LHX9",
    "EMX1", "EMX2", "GBX1", "GBX2",
    "OTX1", "OTX2", "PAX6", "PAX2", "PAX3", "PAX7",
    "SOX1", "SOX2", "SOX3", "SOX4", "SOX5", "SOX6", "SOX9", "SOX10", "SOX11", "SOX13", "SOX21",
    "TFEB", "TFE3", "MITF", "TFEC",
    "ESR1", "ESR2", "NR3C1", "NR3C2", "AR", "PGR", "PPARA", "PPARD", "PPARG",
    "RARA", "RARB", "RARG", "RXRA", "RXRB", "RXRG",
    "RORA", "RORB", "RORC",
    "PBX1", "PBX2", "PBX3", "PBX4", "MEIS1", "MEIS2", "MEIS3",
    "HEY1", "HEY2", "HEYL", "HES1", "HES5", "HES6",
    "ZBTB7A", "ZBTB7B", "ZBTB16", "ZBTB18", "ZBTB20",
    "ASCL1", "ASCL2", "ASCL3", "ASCL4",
    "NGN1", "NGN2", "NGN3", "NEUROG1", "NEUROG2", "NEUROG3",
}

# Order matters: first hit wins. Mechanotransduction is resolved BEFORE
# calcium_ion_channel so dual-membership genes (TRPV4, TRPM7, TRPC6 ...)
# land in the mechanical class. CCN1/CYR61, CCN2/CTGF, JNK/p38 members
# fall to immediate_early_stress / transcription_factor respectively
# because those classes are listed earlier or after as appropriate.
GENE_CLASS_SETS: dict[str, set[str]] = {
    "immediate_early_stress": IMMEDIATE_EARLY_STRESS,
    "mechanotransduction": MECHANOTRANSDUCTION,
    "ecm": ECM,
    "synaptic_plasticity": SYNAPTIC_PLASTICITY,
    "calcium_ion_channel": CALCIUM_ION_CHANNEL,
    "TDP43_RNA_splicing": TDP43_RNA_SPLICING,
    "mitochondrial_nuclear_encoded": MITO_NUCLEAR_LITERALS,  # plus prefix scan below
    "proteostasis": PROTEOSTASIS,
    "apoptosis_stress": APOPTOSIS_STRESS,
    "transcription_factor": TRANSCRIPTION_FACTOR,
}


def classify_gene(symbol: str) -> str:
    """Return the biological gene class for `symbol`, or 'other' if no match."""
    if not symbol:
        return "other"
    for cls, members in GENE_CLASS_SETS.items():
        if symbol in members:
            return cls
    # Prefix scan for nuclear-encoded mitochondrial
    if symbol.startswith(MITO_NUCLEAR_PREFIXES):
        return "mitochondrial_nuclear_encoded"
    return "other"


def is_artifact_symbol(symbol: str) -> tuple[bool, str]:
    """Return (is_artifact, reason). `reason` is empty if not an artifact."""
    if not symbol:
        return False, ""
    if symbol.startswith("MT-"):
        return True, "mitochondrial-encoded"
    if symbol.startswith("RPL") or symbol.startswith("RPS"):
        return True, "cytoplasmic ribosomal protein"
    if symbol in ARTIFACT_LITERALS:
        return True, "abundant lncRNA"
    if symbol in HEMOGLOBINS:
        return True, "hemoglobin"
    if symbol in NONNEURONAL_MARKERS:
        return True, "non-neuronal contamination marker"
    return False, ""
