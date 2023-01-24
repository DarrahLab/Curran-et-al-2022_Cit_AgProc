# Citrullination modulates antigen processing and presentation by revealing cryptic epitopes in rheumatoid arthritis

## Abstract
Cryptic peptides, hidden from the immune system under physiologic conditions, are revealed by changes to MHC class II processing and hypothesized to drive the loss of immune tolerance to self-antigens in autoimmunity. Rheumatoid arthritis (RA) is an autoimmune disease characterized by immune responses to citrullinated self-antigens, in which arginine residues are converted to citrullines. We investigated the hypothesis that citrullination exposes cryptic epitopes by modifying protein structure and proteolytic cleavage. We show that citrullination alters processing and presentation of autoantigens, resulting in the generation of a unique citrullination-dependent repertoire composed primarily of native sequences. This repertoire stimulates T cells from RA patients with anti-citrullinated protein antibodies more robustly than controls. The generation of this unique repertoire is achieved through altered protease cleavage and protein destabilization, rather than direct presentation of citrulline-containing epitopes, suggesting a novel paradigm for the role of protein citrullination in the breach of immune tolerance in RA.

## File Structure 
All raw data necessary to execute our scripts are included in the SourceData subfolder. <br />
  **01_SeqSub** Identifies citrullination sites in MS data and performs substitutions at those sites to prepare input to alphafold. <br />
  **02_DiffExbyRes** Performs differential expression calculations between native, pad2-citrullinated, and pad4-citrullinated mass spec data by residue. <br />
  **03_region.cit.dist** Identifies Created or Destroyed regions of protein antigens and calculates the distance from these regions to the nearest citrullinated residue. <br />
  **04_stats** Calculates statistics on regional distance data from script 03 above. <br />
  **05_plotting** Plots differential expression in the format used in the final manuscript. <br />
  **06_3DStructureAnalysis** Performs distance calculations between Created/Destroyed regions and citrullinated residues on alphafold predictions of structure. <br />
  **07_VolcanoPlots_NetMHCII-cores.rmd** Creates volcano plots and calculates abundance of NetMHCII-predicted high-binding cores in total peptide repertoire. <br />
