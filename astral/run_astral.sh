# for astral 
cd ~/Documents/bom_phy/analyses/astral_multilocus

# make input file 
cat gt_realBS/gt_output/*.tre > in.trees
# OR see the other script for adding new lines 

# OPTIONAL collapse low support branches in gene trees (need newick utils)
nw_ed  in.trees 'i & b<=10' o > in_BS10.trees


# run astral (on Carrie's computer)

# just get best ASTRAL Tree 
java -jar /Applications/Phylogeny_Programs/Astral/astral.5.7.7.jar -i in_BS10.trees -o output/Bomarea_BS10_2.tre 2> output/Bomarea_BS10_2.log
# also do multilocus bootstrapping 
java -jar /Applications/Phylogeny_Programs/Astral/astral.5.7.7.jar -i in_BS10.trees -b bootstraps.txt -o output/Bomarea_multilocus2_BS10.tre 2> output/Bomarea_multilocus2_BS10.log
