# extract the entrez ids from the gtf file (downloaded from NCBI)
gffpath="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0"
gfffile="/mnt/sda1/RNA/40-815970407/Sheep/Reference_geneome/ARS-UI_Ramb_v3.0/genomic.gff"

zgrep 'GeneID' $gfffile \
  | cut -f9 | perl -pe 's/ID.*(GeneID:\d+).*gene=([^;]*).*/\1\t\2/g' \
  | sort -u > $gffpath/extracted_entrezids

# remove GeneID from all the lines and add column names
sed 's/GeneID://g' $gffpath/extracted_entrezids > $gffpath/entrezids 
sed  -i '1i GeneID\tGeneName' $gffpath/entrezids

# getting description
grep 'GeneID' $gfffile \
  | cut -f9 | perl -pe 's/ID.*(GeneID:\d+).*description=([^;]*).*/\1\t\2/g' \
  | sort -u > $gffpath/extracted_entrezids_with_description
  
# remove GeneID from all the lines and add column names
sed 's/GeneID://g' $gffpath/extracted_entrezids_with_description > $gffpath/entrezids_with_description
sed  -i '1i GeneID\tdescription' $gffpath/entrezids_with_description

rm extracted_entrezids_with_description
rm extracted_entrezids

# extract the gene ontology ids and gene ids from the .gaf file (downloaded from NCBI)
awk -F '\t| |"' -v OFS="\t" '{print $2,$5}' $gffpath/GCF_016772045.2-RS_2023_10_gene_ontology.gaf > geneids_GOIds

# combine the GOids and description file
# awk 'NR==FNR {h[$1] = $2; next} {print $1,$2,$3,h[$1]}' $gffpath/entrezids_with_description $gffpath/geneids_GOIds > $gffpath/entrezids_with_GOids_and_description
