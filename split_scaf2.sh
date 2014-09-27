#split_scaf.sh takes a .scafSeq file as input and creates a folder with each individual scaffold in its own .fasta file

#echo "still working0"
#while read LINE;do
#$infile=LINE
#done

#get beginning of file name
file=`awk 'END{print FILENAME}' $1 | awk 'BEGIN{FS="."} {print $1}'`

sed 's/[, |]/_/g' $1 | sed 's/__/_/g' | sed 's/_$//g'  | awk '$1 ~ /^>/ {print $0 "QQQQ"} $1 !~ /^>/ {print $0}'  | tr -d "\n\t ^M" | sed 's/>/\n>/g'  | sed 's/QQQQ/\n/g' | sed '/^$/d' | sed 's/ /_/g' | sed 's/[.]//g' > ${file}_unwrapped.fasta


mkdir split_scaf
grep '>' ${file}_unwrapped.fasta > headerList.sp                                                        #get list of fasta headers
sed 's/>//g' headerList.sp > headerList2.sp             #replace spaces and carats with underscores


headerVar=`cat headerList2.sp`                                                          #read as local variable
headerArray=($headerVar)
numHeaders=`echo ${#headerArray[@]}`

#get beginning of file name
#file=`awk 'END{print FILENAME}' $1 | awk 'BEGIN{FS="."} {print $1}'`

#get rid of all white spaces, except carriage returns between headers and sequences
#sed 's/[, |]/_/g' $1 | sed 's/__/_/g' | sed 's/_$//g'  | awk '$1 ~ /^>/ {print $0 "QQQQ"} $1 !~ /^>/ {print $0}'  | tr -d "\n\t ^M" | sed 's/>/\n>/g'  | sed 's/QQQQ/\n/g' | sed '/^$/d' > ${file}_unwrapped.fasta


for ((i=0; i<$numHeaders; i++))
do
        cHeader=${headerArray[$i]}
        grep -A1 $cHeader ${file}_unwrapped.fasta > split_scaf/${headerArray[$i]} #cast header list as array
done

rm ${file}_unwrapped.fasta
rm *sp
