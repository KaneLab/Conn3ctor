#format fasta file for downstream processes
#$1 in fasta file to correct

#get beginning of file name
file=`awk 'END{print FILENAME}' $1 | awk 'BEGIN{FS="."} {print $1}'`

#edit genbank name formatting, which gums up the rest of the workflow
#get rid of all white spaces, except carriage returns between headers and sequences
sed 's/[, |]/_/g' $1 | sed 's/__/_/g' | sed 's/_$//g'  | awk '$1 ~ /^>/ {print $0 "QQQQ"} $1 !~ /^>/ {print $0}'  | tr -d "\n\t ^M" | sed 's/>/\n>/g'  | sed 's/QQQQ/\n/g' | sed '/^$/d'
