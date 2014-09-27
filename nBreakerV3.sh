sed 's/\(.*\)\s.*/\1/g' $1 > pre_infile #simplify contig names
sed 's/caffold//g' pre_infile > infile
namesO=`grep '>' infile` #get headers
nameArrayO=($namesO)  # make name an array (may not be necessary)
posLengthO=${#nameArrayO[@]}


#get rid of all carriage returns
sed -i ':a;N;$!ba;s/\n//g' infile
for ((i=0; i<$posLengthO; i++)); do
    namy=${nameArrayO[$i]}
    namyPos=`cat -n infile | grep $namy | awk '{print $1}'`
    if [ $namyPos -eq 1 ]; then
        sed -i "s/\($namy\)/\1\n/" infile
    else
        sed -i "s/\($namy\)/\n\1\n/" infile
fi
done

names=`grep -B1 'N' infile | grep '>'`
nameArray=($names)
posLength=${#nameArray[@]}

keyword=">BAD"


sed "s/N\([ATGC]\)/N\n$keyword\n\1/g" infile > outfile #split all scaffolds with N's into multiple scaffolds

replaceN=`cat -n outfile | grep $keyword | awk '{print $1}'` #get line numbers that contain key word
endPos=`cat -n outfile | awk 'END{print $1}'`
replaceArray=() #cast as array
replaceArray+=($replaceN)
replaceArray+=($endPos)


headPosArray=()
headPosArrayArray=()
for (( k=0; k<$posLength; k++ ));do
        currentHead=${nameArray[$k]}
        headPos=`cat -n outfile | grep $currentHead | awk '{print $1}'`
        headPosArrayArray=($headPos)
        headPosArray+=(${headPosArrayArray[0]})
done
headPosArray+=($endPos)


for (( j=0; j<$posLength; j++ ));do

    newHead=${nameArray[$j]}
    newHeadPos=${headPosArray[$j]}

        if [ $newHeadPos -lt ${headPosArray[$(($posLength-1))]} ];then
                nextHeadPos=${headPosArray[$(( $j + 1 ))]}
        else
                endPos=`cat -n outfile | awk 'END{print $1}'`
                nextHeadPos=($endPos)
        fi

    nextJ=$(( $j + 1 ))


        for (( i=0; i<${headPosArray[$nextJ]}; i++));do
            if [ ${replaceArray[$i]} -lt ${headPosArray[$nextJ]} ];then
                replace=${replaceArray[$i]} #make notation easier
                sed -i -e "${replace}s/$keyword/${newHead}_${i}/" outfile
            else
                break
            fi
        done
done
sed -i 's/N//g' outfile
sed -i 's/[ ][ACTG]//g' outfile
rm infile; rm pre_infile
echo "Finished breaking Ns."
