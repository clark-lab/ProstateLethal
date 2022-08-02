#This script used to run GAT software (Heger 2013 Bioinformatics: doi.org/10.1093/bioinformatics/btt343) to find overlaps between DMRs and ChromHMM 

#script uses ChromHMM data generated as part of Giles 2019 Epigenetics and Chromatin: doi.org/10.1186/s13072-019-0258-9  ('metadata/PrEC_15_segments.bed' and 'metadata/LNCaP_15_segments.bed')

#script uses segment and annotation data generated through scripts 03_02_01_MakeDMRBedfiles_forGAT.R and 03_02_02_MakeWGBSBackgroundBedfile_forGAT.R

#script run from the same BASEDIR as used in all other scripts



###### Hyper DMRs

#GAT on PrEC ChromHMM - 15 state
gat-run.py --segments="WGBS/results/hyperDMR_analysis.bed" --annotations="metadata/PrEC_15_segments.bed" --workspace="WGBS/results/AllCpGs_analysis.bed" --nbuckets=10000 --bucket-size=10650 --num-samples=10 > "WGBS/results/GAT_HyperDMR_PrEC15.tsv"

#GAT on LNCaP ChromHMM - 15 state
gat-run.py --segments="WGBS/results/hyperDMR_analysis.bed" --annotations="metadata/LNCaP_15_segments.bed" --workspace="WGBS/results/AllCpGs_analysis.bed" --nbuckets=10000 --bucket-size=10650 --num-samples=10 > "WGBS/results/GAT_HyperDMR_LNCaP15.tsv"



######Hypo DMRs

#GAT on PrEC ChromHMM - 15 state
gat-run.py --segments="WGBS/results/hypoDMR_analysis.bed" --annotations="metadata/PrEC_15_segments.bed" --workspace="WGBS/results/AllCpGs_analysis.bed" --nbuckets=10000 --bucket-size=10650 --num-samples=10 > "WGBS/results/GAT_HypoDMR_PrEC15.tsv"

#GAT on LNCaP ChromHMM - 15 state
gat-run.py --segments="WGBS/results/hypoDMR_analysis.bed" --annotations="metadata/LNCaP_15_segments.bed" --workspace="WGBS/results/AllCpGs_analysis.bed" --nbuckets=10000 --bucket-size=10650 --num-samples=10 > "/WGBS/results/GAT_HypoDMR_LNCaP15.tsv"

 