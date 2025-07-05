
# Load structures
load human_tp53.pdb
load pig_tp53.pdb

# Align structures
align pig_tp53, human_tp53

# Visualize
spectrum b, rainbow, human_tp53
show cartoon, human_tp53
show surface, pig_tp53, transparency=0.6

# Highlight CTD
select human_ctd, resi 360-393
select pig_ctd, resi 360-393
show sticks, human_ctd
show mesh, pig_ctd
color red, human_ctd
color blue, pig_ctd

# Save image
ray 1600, 1200
png structural_comparison.png
