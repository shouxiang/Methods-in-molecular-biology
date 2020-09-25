This repository contains the source codes for flow cytometric analyses in the book chapter titled "Measuring cysteine exposure in unfolded proteins with tetraphenylethene maleimide and its analogues" for an edition on  “The unfolded Protein Response – Methods and Protocols” in the lab protocol series Methods in Molecular Biology. Souce codes are written based on .csv files in the "Examples" folder. Below is the instruction on how to use the codes.

1. .fcs files are converted to .csv files by FlowJo. Place all .csv files in a folder. The "Examples" folder includes published data and can be used as a demo.
2. Run “gating.R” to modify and determine the gating strategy.
3. Modify the gates in “initialize.R”, and run “initialize.R” to pre-process raw .csv files.
4. Run “statistical analysis.R” to compare control and experimental groups and generate plots.
