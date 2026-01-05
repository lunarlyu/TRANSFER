"""
Constants Module

This module contains all configuration constants used across the
cell surface marker prediction pipeline.
"""

# Immune cell subtypes that should be grouped into "immune" tissue
IMMUNE_SUBTYPES = [
    "t-cells", "b-cells", "nk-cells", "monocytes", "macrophages",
    "granulocytes", "plasma cells", "dendritic cells"
]

# Non-membrane genes to exclude (extracted from UniProt data)
# These genes were identified as not being cell surface proteins
NON_MEMBRANE_GENES = {
    'YBX1', 'PPT1', 'C6orf120', 'RPL35', 'AEBP1', 'ANGPT1', 'PRG2', 'LAMB2',
    'LMAN2L', 'NAGPA', 'PLD3', 'IFI30', 'MAN2A2', 'NCLN', 'TOR1AIP2', 'SEMA4C',
    'IL17D', 'LAMB1', 'STS', 'RPL28', 'MAN2B1', 'ATP6AP1', 'CRELD2', 'BRI3BP',
    'ADAM15', 'ORM1', 'GBA', 'LIPA', 'SCARA3', 'SLC5A3', 'LTBP1', 'PXDN',
    'ST8SIA4', 'B3GNT2', 'HEPACAM2', 'NUP210', 'CTSC', 'PFN1', 'GALNS',
    'TUBB2C', 'KNG1', 'MOXD1', 'WNT5B', 'FKTN', 'WNT11', 'ELOVL4', 'POSTN',
    'GDF15', 'SNRPD2', 'WNT5A', 'SLC43A3', 'ABI3BP', 'GRN', 'EPX', 'SERPINI1',
    'PLD4', 'APOB', 'LMNA', 'FAM20B', 'CRB2', 'RPL35A', 'UBA1', 'LZTFL1',
    'LAMA5', 'RPL9', 'LTF', 'EMILIN1', 'COL14A1', 'PPT2', 'MIA3', 'ATRX',
    'HGSNAT', 'COL11A1', 'HSP90B1', 'APOH', 'BTD', 'TFPI2', 'VWF', 'NXPH4',
    'PLAU', 'FBN2', 'PCBP2', 'CFH', 'FKBP9', 'METTL9', 'HNRNPU', 'ATF6B',
    'NPTX1', 'C4BPB', 'SV2A', 'TUBB', 'PPIB', 'ITPR3', 'ANGPTL4', 'CHGA',
    'CD46', 'OGFOD1', 'DCN', 'GDPD5', 'LAMC2', 'F5', 'PRELP', 'ASAH1',
    'QSOX1', 'ISLR', 'GMPS', 'TGFB1', 'PLOD3', 'GMPPB', 'TMEM9', 'GLB1',
    'APOD', 'ITIH2', 'GAPDHS', 'LGMN', 'IGLON5', 'CASD1', 'SERPINB2',
    'SERPINF1', 'NXF2', 'OGN', 'DMXL2', 'FN1', 'BDNF', 'PRSS23', 'PLBD1',
    'RCN3', 'RPL18', 'CCDC80', 'ELANE', 'HEXA', 'ABCC5', 'PON1', 'HP',
    'IGF2R', 'SULF2', 'RPL21', 'TIMP1', 'YWHAQ', 'CCDC126', 'TOR1B',
    'PTGFRN', 'TPP1', 'COL6A3', 'HSPA13', 'ABCC4', 'ACTN4', 'RNF13', 'FGL2',
    'COL6A1', 'PCBP3', 'DKK3', 'FGB', 'BCHE', 'RSL1D1', 'NAGA', 'LEMD2',
    'CHSY3', 'SERPINA1', 'GALNT5', 'ADAMTS15', 'POFUT2', 'HSD11B1', 'COL3A1',
    'TOR3A', 'HMCN1', 'C4BPA', 'CALR', 'NOP56', 'GNPTAB', 'SERPINH1',
    'CASP7', 'THSD4', 'PLAT', 'VCAN', 'EDEM3', 'CGN', 'RPL17', 'PSAP',
    'MAN2B2', 'LMF2', 'TLL2', 'TLL1', 'GNS', 'CTSD', 'ERAP1', 'TYR', 'HYOU1',
    'PTPRQ', 'RPL3', 'HEXB', 'F2', 'ARSF', 'GNPTG', 'TF', 'COL12A1',
    'SLCO3A1', 'STT3B', 'OLFML2B', 'COL1A1', 'EXTL2', 'GZMH', 'FKBP14',
    'PCYOX1', 'M6PR', 'SEMA3C', 'TMX4', 'ARSB', 'COL6A6', 'LMAN2', 'PPIAL4A',
    'LRPAP1', 'DNASE1L1', 'HSPD1', 'SPANXB1', 'INTS12', 'APLF', 'TUBA4A',
    'THBS1', 'TNFRSF11B', 'B4GALT3', 'EFEMP1', 'CDSN', 'COL5A2', 'RPL32',
    'CHST3', 'LTBP3', 'WDR13', 'LAMA1', 'HNRNPK', 'FKBP7', 'MPO', 'MMP13',
    'ERGIC2', 'ANGPTL1', 'ACP2', 'SLC29A2', 'GLCE', 'CTSA', 'ADAMTSL4',
    'OLFML3', 'UST', 'TRAP1', 'CERCAM', 'LGALS3BP', 'SERPING1', 'NPTX2',
    'FGA', 'CTSL1', 'INHBE', 'IGFBP3', 'UGGT1', 'PBXIP1', 'B4GALT5', 'SEL1L',
    'GLA', 'HSD17B2', 'MAN2A1', 'ANGPTL2', 'SULF1', 'FLNA', 'LAMC1', 'CES1',
    'BIRC6', 'SUMF1', 'ACTB', 'FLG2', 'SLIT2', 'SLC22A4', 'RPL37A', 'SLIT3',
    'RPS23', 'LOXL4', 'PIGT', 'GPI', 'ENTPD7', 'SERINC1', 'LUM', 'OSTM1',
    'OLFML2A', 'LIPG', 'A2M', 'COL2A1', 'FBLN2', 'FUT11', 'TPI1', 'HPX',
    'C4A', 'SEMA3A', 'SFRS3', 'CCT4', 'SERPINE2', 'CHST7', 'TXNDC15', 'RALY',
    'LAMP3', 'PI16', 'FICD', 'SLCO4C1', 'UQCRC2', 'CREG1', 'EEF1D', 'RPL27A',
    'FSTL1', 'CRTAP', 'DPP7', 'NAV2', 'ZNF844', 'VEGFC', 'RCN1', 'TFPI',
    'DMTF1', 'ASPN', 'P4HA2', 'MPHOSPH9', 'XRCC6', 'NELL1', 'ASPH', 'NTN3',
    'TNC', 'RPS15A', 'AGPS', 'MYH9', 'PRSS35', 'COL5A1', 'ABCC2', 'SERPINA3',
    'VTN', 'TUBB6', 'DBH', 'SMPD1', 'CPZ', 'TMED9', 'CP', 'MXRA5', 'AKAP11',
    'MMRN1', 'ERAP2', 'COL7A1', 'NAGLU', 'PLTP', 'MFAP4', 'MGAT4B', 'PRKCSH',
    'C1QTNF2', 'ZC3HAV1', 'P4HA1', 'TNXB', 'EMILIN2', 'OLFM1', 'SIL1',
    'CRLF1', 'SCARB2', 'GREM1', 'PA2G4', 'RPL19', 'UGGT2', 'FKBP10',
    'COL5A3', 'IDH2', 'EPDR1', 'PLOD2', 'BGN', 'RERG', 'PLBD2', 'FBN1',
    'PTBP1', 'AHSG', 'ACAA1', 'PLOD1', 'GUSB', 'BMP1', 'LOX', 'MFAP5',
    'LAMB3', 'C3', 'ST6GALNAC3', 'COL1A2', 'F13A1'
}

# Known positive markers from literature (for validation during optimization)
# Format: [gene_name, tissue, cell_type, delivery_method]
LITERATURE_POSITIVES = [
    ["PECAM1", "lung", "endothelial cells", "LNP"],
    ["VCAM1", "vascular", "endothelial cells", "LNP"],
    ["CD4", "pbmc", "t-cells", "LNP"],
    ["CD5", "pbmc", "t-cells", "LNP"],
    ["CD19", "pbmc", "b-cells", "LNP"],
    ["CD3", "pbmc", "t-cells", "LNP"],
    ["NCR1", "pbmc", "nk-cells", "LNP"],
    ["CD14", "pbmc", "macrophages", "LNP"],
    ["MRC1", "pbmc", "macrophages", "LNP"],
    ["ITGAM", "pbmc", "macrophages", ""],
    ["CD28", "pbmc", "t-cells", "EDV"],
    ["CD40", "pbmc", "b-cells", "lenti"],
    ["ENG", "vascular", "endothelial cells", "LNP"],
    ["MRC1", "pbmc", "dendritic cells", "LNP"],
    ["CD8", "pbmc", "t-cells", "LNP"],
    ["PDPN", "skin", "endothelial cells", "LNP"],
    ["PLVAP", "lung", "endothelial cells", "LNP?"],
    ["FCER2", "pbmc", "b-cells", ""],
]

# LITERATURE_POSITIVES with immune cells grouped into "immune" tissue
# This matches the tissue naming convention used in the processed data
# Generated automatically by run_data_processing() - do not edit manually
# LITERATURE_POSITIVES with immune cells grouped into "immune" tissue
# This matches the tissue naming convention used in the processed data
# Generated automatically by run_data_processing()
LITERATURE_POSITIVES_GROUPED = [
    ['PECAM1', 'lung', 'endothelial cells', 'LNP'],
    ['VCAM1', 'vascular', 'endothelial cells', 'LNP'],
    ['CD4', 'immune', 't-cells', 'LNP'],
    ['CD5', 'immune', 't-cells', 'LNP'],
    ['CD19', 'immune', 'b-cells', 'LNP'],
    ['CD3', 'immune', 't-cells', 'LNP'],
    ['NCR1', 'immune', 'nk-cells', 'LNP'],
    ['CD14', 'immune', 'macrophages', 'LNP'],
    ['MRC1', 'immune', 'macrophages', 'LNP'],
    ['ITGAM', 'immune', 'macrophages', ''],
    ['CD28', 'immune', 't-cells', 'EDV'],
    ['CD40', 'immune', 'b-cells', 'lenti'],
    ['ENG', 'vascular', 'endothelial cells', 'LNP'],
    ['MRC1', 'immune', 'dendritic cells', 'LNP'],
    ['CD8', 'immune', 't-cells', 'LNP'],
    ['PDPN', 'skin', 'endothelial cells', 'LNP'],
    ['PLVAP', 'lung', 'endothelial cells', 'LNP?'],
    ['FCER2', 'immune', 'b-cells', ''],
]

# Genes used as cutoff markers for trimming recommendations
# Markers appearing after these genes in the ranked list are excluded
CUTOFF_GENES = ["SIGLEC5", "PCDHGC5", "PCDHGC11"]

# Cell types to exclude from final recommendations
# These are typically low-quality or redundant cell type annotations
EXCLUDE_CELLS = [
    "adipose tissue fibroblasts",
    "adipose tissue smooth muscle cells",
    "breast adipocytes",
    "breast breast glandular cells",
    "breast breast myoepithelial cells",
    "breast endothelial cells",
    "breast smooth muscle cells",
    "bronchus basal respiratory cells",
    "bronchus club cells",
    "bronchus ionocytes",
    "bronchus smooth muscle cells",
    "colon enteroendocrine cells",
    "colon intestinal goblet cells",
    "colon paneth cells",
    "colon undifferentiated cells",
    "endometrium glandular and luminal cells",
    "endometrium smooth muscle cells",
    "esophagus basal squamous epithelial cells",
    "eye bipolar cells",
]

