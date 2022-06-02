#!/bin/sh

# Split Big Fasta by ID , run Igblast and Get Germline sequece
# Author: DONG Chuang CIML/PMlab 24.07.2020

if [ $# -lt 1 ] ; then
        echo  "\nAuthor: DONG Chuang CIML/PMlab 24.07.2020"
        echo  "USAGE: $0 <command>"
        echo  "For Example     $0 Output/2_FastaByID\n"
        exit 1
fi

mkdir -p $1
cd $1
here=$PWD
ids=( "CLONOTYPE10/CATGCCTAGCCAGTTT" 
      "CLONOTYPE140/CATTATCAGCTGAACG" 
      "CLONOTYPE196/CAGTCCTAGTACACCT"
      "CLONOTYPE234/GCTGCAGAGCTAACAA"
      "CLONOTYPE31/ACGGGTCGTTTGACTG"
      "CLONOTYPE32/GGCTGGTAGAAACCAT"
      "CLONOTYPE34/TGTTCCGAGGAATTAC"
      "CLONOTYPE36/GAGCAGAGTGTTGGGA"
      "CLONOTYPE67/CTAACTTTCCCACTTG"
      "CLONOTYPE83/GTGGGTCCAGATTGCT" )


for id in "${ids[@]}"
	do
		### Check for dir, if not found create it using the mkdir ##
		[ ! -d "$id" ] && mkdir -p "$id"
		cd $id
		uniqCellid=`echo $id | cut -d"/" -f2`
		echo $uniqCellid > $uniqCellid.ids
		fold1="IGH"
		mkdir $fold1;cd $fold1
		eval "bash /mnt/NAS7/PMlab/moFluMemB_LungClones_GermlineGCtree/Script/extractSeqbyID.sh ../$uniqCellid.ids \\
		  /mnt/NAS7/PMlab/moFluMemB_LungClones_GermlineGCtree/Output/1_FastaAll/Root_IGH.fasta IGH.fasta"
		wait
		`singularity exec /mnt/NAS7/Workspace/dongc/changeO/singularity/mglab_moflumemb_changeo.img \\
      /usr/local/bin/AssignGenes.py igblast -s IGH.fasta -b /mnt/NAS7/PMlab/Germline_Mouse/share/igblast --organism mouse --loci ig --format blast`
    sleep 1s
    perl /mnt/NAS7/PMlab/Germline_Mouse/Script/CreateGermline_Mouse.pl IGH_igblast.fmt7 IGH.fasta out.fasta
		cd ..
		fold2="IGKL"
		mkdir $fold2;cd $fold2
		`bash /mnt/NAS7/PMlab/moFluMemB_LungClones_GermlineGCtree/Script/extractSeqbyID.sh ../$uniqCellid.ids \\
		  /mnt/NAS7/PMlab/moFluMemB_LungClones_GermlineGCtree/Output/1_FastaAll/Root_IGKL.fasta IGKL.fasta`
		wait
		`singularity exec /mnt/NAS7/Workspace/dongc/changeO/singularity/mglab_moflumemb_changeo.img \\
      /usr/local/bin/AssignGenes.py igblast -s IGKL.fasta -b /mnt/NAS7/PMlab/Germline_Mouse/share/igblast --organism mouse --loci ig --format blast`
    sleep 1s  
    perl /mnt/NAS7/PMlab/Germline_Mouse/Script/CreateGermline_Mouse.pl IGKL_igblast.fmt7 IGKL.fasta out.fasta  
		cd $here
	done
