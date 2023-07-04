#***************************************#
# FOOTPRINT IDENTIFICATION              #
# script: 14.2_footprint_wellington.sh  #
#***************************************#


# 14.2 Footprint - Wellington (https://pythonhosted.org/pyDNase/tutorial.html#quick-and-easy-footprinting)


#File_Input (BAM)
Input_BAM=$1
#File_Input (BAM)
Input_BED=$2
#Output directory
Out_Dir=$3

#wellington_footprints.py 
/python/bin/wellington_footprints.py $Input_BED  $Input_BAM $Out_Dir
