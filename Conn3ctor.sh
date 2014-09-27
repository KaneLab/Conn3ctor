# Connector attempts to connect contigs and scaffolds by searching the .fasta files containing all of the short-reads for sequences matching a query near the end of the contig, aligning some of the matches, building a consensus, then using the consensus to extend the contig until it overlaps with another contig in the assembly.

#Dependencies:
# split_scaf.sh
# median2.sh
# rc.pl

#Mandatory Inputs:
# -d sets the grep search depth
# -o:  If -o=1, Conn3ctor orders the contigs in your .scafSeq input by size, working from the smallest contig first.  If -o=0, it starts from the largest contig first.  Any other number will disregard order.
# -n sets displacement, which is distance from end of workingContig (or front, if -b is selected) to start tiling
# -t is the number of iterations we run the addSequence before checking for matches elsewhere in the assembly.
# -k sets the length of query sequence to be extended
# -q sets the quality threshold to add on new sequence from the consensus
# -l is the minimum contig size that we're willing to extend.  This is useful for instances in which there are large numbers of small, redundant contigs.

#initialize Mandatory Inputs:
queryLength=`grep 'queryLength' $1 | cut -f2`
disp=`grep 'disp' $1 | cut -f2`
searchDepth=`grep 'searchDepth' $1 | cut -f2`
qualThresh=`grep 'qualThresh' $1 | cut -f2`
ordered=`grep 'ordered' $1 | cut -f2`
checkTime=`grep 'checkTime' $1 | cut -f2`
minContigLength=`grep 'minContigLength' $1 | cut -f2`


#Website with bash getopts tutorial:  http://wiki.bash-hackers.org/howto/getopts_tutorial
#while getopts ":k:n:d:q:t:l:o:r" opt; do
#         case $opt in
#           k)                                                          # -k sets the length of query sequence to be extended
#               if [[ $OPTARG = -* ]]; then
#                       ((OPTIND--))
#                       continue
#               fi
#               queryLength=$OPTARG
#               ;;
#           n)                                                          # -n sets displacement, which is distance from end of workingContig (or front, if -b is selected)
#               if [[ $OPTARG = -* ]]; then
#                       ((OPTIND--))
#                       continue
#               fi
#               disp=$OPTARG
#               ;;
#           d)                                                          # -d sets the grep search depth
#               if [[ $OPTARG = -* ]]; then
#                       ((OPTIND--))
#                       continue
#               fi
#               searchDepth=$OPTARG
#               ;;
#           q)                                                          # -q sets the quality threshold to add on new sequence from the consensus
#               if [[ $OPTARG = -* ]]; then
#                       ((OPTIND--))
#                       continue
#               fi
#               qualThresh=$OPTARG
#               ;;
#           o)                                                          # If -O=1, Conn3ctor orders the contigs in your .scafSeq input by size, working from the smallest contig first.  If -O=0, it starts from the largest contig first.  Any other number will disregard order.
#               if [[ $OPTARG = -* ]]; then
#                       ((OPTIND--))
#                       continue
#               fi
#               ordered=$OPTARG
#               ;;
#           t)                                                          # -t is the number of iterations we run the addSequence before checking for matches elsewhere in the assembly.
#               if [[ $OPTARG = -* ]]; then
#                       ((OPTIND--))
#                       continue
#               fi
#               checkTime=$OPTARG
#
#               ;;
#           l)                                                          # -l is the minimum contig size that we're willing to extend.  This is useful for instances in which there are large numbers of small, redundant contigs.
#               if [[ $OPTARG = -* ]]; then
#                       ((OPTIND--))
#                       continue
#               fi
#               minContigLength=$OPTARG
#               ;;
#           r)                                                          # -r flag removes all files in /connected, /extended/, failed/, and split_scaf/ .
#               remove=$OPTARG
#               ;;
#          \?)
#               echo "Invalid Argument.  Variables must be -d -o -n -t -k -q -l"
#               exit
#          ;;
#        esac
#done

echo "Type 'y' to remove connected/*, extended/*, failed/*, split_scaf/*, quarantine/*, smallContigs/*, and temp_connected/*."
read removeFlag

if [ $removeFlag == "y" ];then
        rm -r connected/ ; rm -r extended/ ; rm -r failed/ ; rm -r quarantine/ ; rm -r split_scaf/; rm -r smallContigs/; rm -r temp_connected
elif [ $removeFlag == "yes" ];then
        echo "I told you to just use 'y', dum-dum!  Okay, just this once."
        rm -r connected/ ; rm -r extended/ ; rm -r failed/ ; rm -r quarantine/ ; rm -r split_scaf/; rm -r smallContigs/; rm -r temp_connected
fi


mkdir quarantine; mkdir failed; mkdir smallContigs; mkdir connected; mkdir extended; mkdir temp_connected
echo " " > connected/empty
echo " " > extended/empty
echo " " > failed/empty
echo " " > temp_connected/empty
echo " " > log

#echo "Dollar 15 is ${15}"
#bash split_scaf2.sh ${15} #$15 means the first positional argument, which is the infile. Infile MUST BE AT 15th POSITION!!!
numArg=$#
bash split_scaf2.sh ${!numArg}

#call split_scaf.sh
if [ $ordered -eq 1 ];then
        ls split_scaf -Sor > scafList
        awk -v s=$minContigLength '$4 > s {print $8}' scafList > scafList.co
elif [ $ordered -eq 0 ];then
        ls split_scaf -So > scafList
        awk -v s=$minContigLength '$4 > s {print $8}' scafList > scafList.co
else
        ls split_scaf > scafList
        awk -v s=$minContigLength '$4 > s {print $8}' scafList > scafList.co
fi




cat scafList.co
scafString=`cat scafList.co`
scafListArray=($scafString)
scafListLength=${#scafListArray[@]}
firstLoop=1
universalAddCount=0

cat scafList.co
scafString=`ls split_scaf/`
scafListArrayL=($scafString)
scafListLengthL=${#scafListArrayL[@]}
firstLoop=1
universalAddCount=0



#format files in split_scaf to only have two lines: a header, and all of the sequence on the next line (with carriage returns removed)

for ((h=0; h<$scafListLengthL; h++));do
        header1=${scafListArrayL[$h]}
        header1Size=`ls -So split_scaf/$header1 | awk '{print $4}'`
        #header1SizeT2=($header1SizeT)
        #echo "header1size2 is $header1SizeT2"
        #header1Size=${header1SizeT2[0]}
        #echo "header size is ${header1SizeT2[0]}"
        if [ $header1Size -ge $minContigLength ]; then
                tempString=`grep -v '>' ./split_scaf/$header1 | sed ':a;N;$!ba;s/\r//g'| sed ':a;N;$!ba;s/\n//g' | sed 's/ //g'`
                rm ./split_scaf/$header1
                echo '>'$header1 > ./split_scaf/$header1
                echo $tempString >> ./split_scaf/$header1
        else
                mv split_scaf/$header1 smallContigs/
        fi
done


#MEGA LOOOOOP!!!!
for ((i=0; i<$scafListLength; i++));do

checkCounter=$checkTime

#Format the input contig or scaffold to remove spaces, and carriage returns, then return it back into fasta format by adding back in the header.
        header1=${scafListArray[$i]}
        grep -v '>' ./split_scaf/$header1 | sed ':a;N;$!ba;s/\n//g' | sed 's/ //g' > workingContig.co
        echo "Moved $header1 to quarantine"
        mv split_scaf/$header1 quarantine/
#Get the working contig, its length.
        workingContig=`grep -v '>' workingContig.co`
        originalContigLength=`echo $((${#workingContig} - $disp))`


globalAddCount=0


while [ $checkCounter -gt 0 ];do
        if [ $firstLoop -eq 0 ];then
        contigLength=`echo ${#workingContig}`
        else
        preContigLength=`echo ${#workingContig}`
        contigLength=$(( $preContigLength - $disp ))
        fi


######################Start addSequence preparation.  This includes getQuery, getAlignment, checkNumConcensus, formatConcensus, stringToArray, and truncateWorkingContig##############################


#getQuery:  This is based on values of -k (queryLength) and -n (disp).
        if [ $firstLoop -eq 0 ];then
                query=`echo ${workingContig:$contigLength-$queryLength:$queryLength}`
        else
                query=`echo ${workingContig:$contigLength-$disp-$queryLength:$queryLength}`
        fi
        firstLoop=0

#getAlignment:  grep the query in the short-reads and save them.  Note that $searchDepth is the number of matching lines that grep takes FROM EACH FILE.  So max number of reads = $searchDepth*(number of files in directory).
        grep -m $searchDepth -h $query ./reads/*fastq > reads.co
        numReads=`wc -l reads.co | awk '{ print $1 }'`
        if [ $numReads -eq 0 ];then             #if your grep query came up empty, move on to next contig
                echo "No reads matched your grep query!  Dammit, Jim!  I'm a conn3ctor, not a miracle worker!"
                echo '>'$header1 >> failed/$header1
                echo $workingContig >> failed/$header1
                rm ./quarantine/$header1
                firstLoop=1
                checkCounter=0
                break
        fi
        cat reads.co | awk '{ print ">" FNR "\n" $0 }' > reads_formatted.co             #format output to look like fasta
        cap3 reads_formatted.co                                                 #perform cap3 Alignment on reads


#checkNumConcensus:  If cap3 finds more than one consensus, we need an adult.
        numConsensus=`grep -c '>' reads_formatted.co.cap.contigs`
        if [ $numConsensus -gt 1 ]; then
                median=
                lineNumstring=`cat -n reads_formatted.co.cap.contigs.qual | grep '>' | awk '{print $1}'`
                lineNums=($lineNumstring)
                totalLinestring=`wc -l reads_formatted.co.cap.contigs.qual`
                totalLines=($totalLinestring)
                for ((k=0; k<=$(($numConsensus-1)); k++));do
                if [ $k -eq $(($numConsensus-1)) ];then         #Need special condition for the last consensus in the file...
                        line1=${lineNums[$k]}
                        cat reads_formatted.co.cap.contigs.qual | sed -n "$line1,${totalLines}p" > consensus.co
                        median+=`bash median2.sh consensus.co | awk '{print " " $1}'`
                        medianArray=($median)
                else                                                    #Find median of each consensus.
                        line1=${lineNums[$k]}
                        line2=$((${lineNums[$(($k+1))]} - 1 ))
                        cat reads_formatted.co.cap.contigs.qual | sed -n "$line1,${line2}p" > consensus.co

                        median+=`bash median2.sh consensus.co | awk '{print " " $1}'`
                        medianArray=($median)
                        fi
                done
                echo "More than one consensus sequence.  Choosing contig with highest median."
        max=0
        arrayLength=`echo ${#medianArray[@]}`
                for ((n=0; n<arrayLength; n++));do                      #Find the array with the highest median.
                        if [ ${medianArray[$n]} -ge $max ];then
                                max=${medianArray[$n]}
                                maxPos=$n
                        fi
                done


                contigLineNumstring=`cat -n reads_formatted.co.cap.contigs | grep '>' | awk '{print $1}'`
                contigLineNums=($contigLineNumstring)
                contigTotalLinestring=`wc -l reads_formatted.co.cap.contigs`
                contigTotalLines=($contigTotalLinestring)
                start=${lineNums[$(($maxPos))]}
                contigStart=${contigLineNums[$(($maxPos))]}
                if [ $(($maxPos+1)) -eq $arrayLength ];then
                        end=$totalLines
                        contigEnd=$contigTotalLines
                else
                        maxPosPlus=$(($maxPos + 1))
                        end=$((${lineNums[$maxPosPlus]} - 1))
                        preContigEnd=${contigLineNums[$maxPosPlus]}
                        contigEnd=$(($preContigEnd - 1))
                fi
                cat reads_formatted.co.cap.contigs.qual | sed -n "$start,${end}p" > quals.co    #Get quals from correct consensus
                cat reads_formatted.co.cap.contigs | sed -n "$contigStart,${contigEnd}p" > correct_consensus.co         #Get correct consensus
                grep -v '>' correct_consensus.co | sed ':a;N;$!ba;s/\n//g'> reads_consensus.co
                cat quals.co    #display correct consensus if several are found

        #FAILED FOLDER:  If no consensus found in cap3 alignment
        elif [ $numConsensus -eq 0 ]; then
                echo "NO CONSENSUS FOUND.  UNABLE TO PROCEED.  PROTOCOL 6.022b CALLS FOR TERMINATION.  HAVE A GOOD DAY."
                echo "$header1:     Moved to /failed for failure to make consensus with query $query." >> log
                #echo "Original contig is $workingContig" >> failed/$header1
                echo '>'$header1 >> failed/$header1
                echo $newSeq >> failed/$header1
                rm ./quarantine/$header1
                firstLoop=1
                checkCounter=0
                break
        else
                maxPos=0
                grep -v '>' reads_formatted.co.cap.contigs.qual > quals.co      #Only one consensus to get quals from
                grep -v '>' reads_formatted.co.cap.contigs | sed ':a;N;$!ba;s/\n//g'> reads_consensus.co        #Only one consensus
                cat quals.co    #display consensus even if only one is found
        fi


#formatConcensus
        qualString=`cat quals.co`                                                       #grab quality scores from cap3 output
        sed -e 's:\(.\):\1 :g' < reads_consensus.co > reads_consensus_spaced.co         #insert space between each nucleotide


#stringToArray:  Consensus sequence and corresponding quality scores turned into vectors
        quals=($qualString)                                                             #string to array
        consensusString=`cat reads_consensus_spaced.co`                                 #grab consensus from cap3 output
        consensus=($consensusString)                                                    #string to array
        consensusLength=`echo ${#consensus[@]}`

        grep -b -o $query reads_consensus.co > position.co                      #find start position of original grep query in the consensus, then print position to file
        queryPos=`cut -d ":" position.co -f1`
        queryPosArray=($queryPos)

        #FAILED FOLDER:  If more than one match to query in cap3 consensus
        if [ ${#queryPosArray[@]} -gt 1 ];then
                if [ `ls failed/ | wc -l` -gt 1 ];then
                        rm failed/empty
                fi
                #echo "The query is $query"
                echo "Multiple matches to query in consensus!  Abort mission!"
                echo "$header1:     Moved to ./failed/ for having more than one match to query $query in consensus." >> log
                #echo "Segmentation faulted with query:   $query" > failed/$header1
                echo '>'$header1 >> failed/$header1
                #echo "Original contig is $workingContig" >> failed/$header1
                echo $newSeq >> failed/$header1
                rm ./quarantine/$header1
                firstLoop=1
                checkCounter=0
                break
        fi
        pos=$(($queryPos + $queryLength))
        addCount=0


#truncateWorkingContig:  Truncate workingContig to only include everything up to the end of the query.  Get rid of sequence in the $disp displacement.
        newSeq=`echo ${workingContig:0:$contigLength}`

heterozygote="false"
#'addSequence' onto $newSeq
        for ((j=$pos; j<=$consensusLength; j++));do
                if [ ${quals[$j]} -ge $qualThresh ];then        #add nucleotides until quality threshold is no longer met
                        newSeq+=${consensus[$pos]}
                        pos=$(( $pos + 1 ))
                        addCount=$(( $addCount + 1 ))
                elif [ $j -lt $consensusLength ] && [ ${quals[$(( $j + 1 ))]} -ge $qualThresh ];then
                        newSeq+=${consensus[$pos]}
                        pos=$(( $pos + 1 ))
                        addCount=$(( $addCount + 1 ))
                        heterozygote="true"
                else
                        checkCounter=$(( $checkCounter - 1 ))           #decrement checkCounter each addSequence loop
                        #EXTENDED FOLDER:  Add contig and any new sequence into ./extended/ if no new sequence can be added and no matches have been found.
                        if [ $addCount -eq 0 ];then
                                echo "The add count for this addSequence run is $addCount.  Nothing else I can do, Dave.  Adding on new sequence and moving on to next contig."
                                echo "$header1:     Added $(( $contigLength - $originalContigLength )) before stalling and heterozygote=$heterozygote."  >> log
                                echo '>'$header1 > extended/$header1
                                echo $newSeq >> extended/$header1
                                #echo "Added $(( $contigLength - $originalContigLength )) nucleotides to $header1." >> ./extended/$header1
                                rm ./quarantine/$header1
                                rm *.co
                                firstLoop=1
                                checkCounter=0
                                #universalAddCount= $(( $universalAddCount + $(( $(( ${#newSeq} - $originalContigLength )) - $disp )) ))
                                break
                        fi
                        echo "Add count is $addCount and check counter is $checkCounter"
                        workingContig=$newSeq
                        if [ $checkCounter -eq 0 ];then         #check every 't' addSeq iterations for conn3ctions
                                echo '>'$header1 > temp_seq.co
                                echo $newSeq >> temp_seq.co
                                matchQuery=${workingContig:${#newSeq}-100:100}
                                matchQueryRC=`perl rc.pl $matchQuery`           #get reverse complement of $matchQuery
                                numMatches=`grep -h -c  $matchQuery ./temp_connected/* ./split_scaf/* ./extended/* ./failed/* | awk '{sum = sum + $1} END {print sum}'`
                                numMatchesRC=`grep -h -c  $matchQueryRC ./temp_connected/* ./split_scaf/* ./extended/* ./failed/* | awk '{sum = sum + $1} END {print sum}'`
                                matchPath=()
                                matchPathRC=()
                                matchSize=()
                                matchPath+=(`grep -l  $matchQuery ./temp_connected/* ./split_scaf/* ./extended/* ./failed/*`)   #find the file path to the matches
                                matchPathRC+=(`grep -l  $matchQueryRC ./temp_connected/* ./split_scaf/* ./extended/* ./failed/*`)
                                biggestMatch=0
                                if [ $(( $numMatches + $numMatchesRC )) -gt 0 ];then
                                        for ((l=0; l<$numMatches; l++));do              #Select the biggest matching contig
                                                matchSize+=(`wc -c ${matchPath[$l]} | awk '{print $1}'`)
                                                if [ ${matchSize[$l]} -gt $biggestMatch ];then
                                                        biggestMatch=${matchSize[$l]}
                                                        correctMatchPath=${matchPath[$l]}
                                                fi
                                        done
                                        for ((m=0; m<$numMatchesRC; m++));do
                                                matchSize+=(`wc -c ${matchPathRC[$m]} | awk '{print $1}'`)
                                                if [ $numMatches -gt 0 ];then
                                                        if [ ${matchSize[$(( $l + $m + 1))]} -gt $biggestMatch ];then
                                                                biggestMatch=${matchSize[$(( $l + $m + 1))]}
                                                                correctMatchPath=${matchPathRC[$m]}
                                                                numMatches=0                    #Allows program to go into the 'elif' loop for joining matches on complementary strand
                                                        fi
                                                elif [ $numMatchesRC -gt 0 ];then
                                                        if [ ${matchSize[$m]} -gt $biggestMatch ];then
                                                                biggestMatch=${matchSize[$m]}
                                                                correctMatchPath=${matchPathRC[$m]}
                                                                numMatches=0                    #Allows program to go into the 'elif' loop for joining matches on complementary strand
                                                        fi
                                                fi
                                        done
                                fi
                                if [ $numMatches -gt 0 ];then   #join conn3ction, if any
                                        match=`cat $correctMatchPath | head -n 1 | xargs -n1 basename | sed 's/>//g'`   #get matching filename
                                        echo "Congratulations!  Added $(( $contigLength - $originalContigLength )) nucleotides to $header1 and joined to contig $match!"
                                        echo "$header1:     Added $(( $contigLength - $originalContigLength )) and joined to $match.  Heterozygote=$heterozygote.  The match was in $correctMatchPath." >> log
                                        echo '>'${header1}_$match >> ./connected/${header1}_$match
                                        echo '>'$header1 > ./temp_connected/$header1
                                        echo $newSeq >> ./temp_connected/$header1
                                        echo $newSeq >> ./connected/${header1}_$match
                                        echo "The query that matches both contigs is: $matchQuery" >> ./connected/${header1}_$match             #include matchQuery to make joining easier
                                        echo "Added $(( $contigLength - $originalContigLength )) nucleotides to $header1 and connected it to $match." >> ./connected/${header1}_$match
                                        cat $correctMatchPath >> ./connected/${header1}_$match
                                        rm ./quarantine/$header1
                                        rm *.co
                                        j=$(( $consensusLength + 1 ))
                                        firstLoop=1
                                        checkCounter=0
                                        #universalAddCount= $(( $universalAddCount + $(( ${#newSeq} - $originalContigLength - $disp )) ))
                                elif [ $numMatchesRC -gt 0 ];then       #Look for conn3tions on reverse strand
                                        matchRC=`cat $correctMatchPath | head -n 1 | xargs -n1 basename | sed 's/>//g'`
                                        echo "Congratulations!  Added $(( $contigLength - $originalContigLength )) nucleotides to $header1 and joined to the reverse complemeneted contig $matchRC!"
                                        echo "$header1:     Added $(( $contigLength - $originalContigLength )) and joined to the reverse complement of $matchRC.  Heterozygote=$heterozygote.  The match was in $correctMatchPath." >> log
                                        echo '>'${header1}_${matchRC}_RC >> ./connected/${header1}_${matchRC}_RC
                                        echo $newSeq >> ./connected/${header1}_${matchRC}_RC
                                        echo '>'$header1 > ./temp_connected/$header1
                                        echo $newSeq >> ./temp_connected/$header1
                                        echo "The query that matches both contigs is: $matchQuery" >> ./connected/${header1}_${matchRC}_RC              #include matchQuery to make joining easier
                                        echo "Added $(( $contigLength - $originalContigLength )) nucleotides to $header1 and connected it to $matchRC." >> ./connected/${header1}_${matchRC}_RC
                                        echo $matchRC >> ./connected/${header1}_${matchRC}_RC
                                        matchPreRC=`tail -n 1 $correctMatchPath`
                                        perl rc.pl $matchPreRC >> ./connected/${header1}_${matchRC}_RC  #Revcomp $matchPreRC and append to file
                                        rm ./quarantine/$header1
                                        rm *.co
                                        j=$(( $consensusLength + 1 ))
                                        firstLoop=1
                                        checkCounter=0
                                        #universalAddCount= $(( $universalAddCount + $(( ${#newSeq} - $originalContigLength - $disp )) ))
                                else checkCounter=$checkTime
                                        echo "No conn3ctions yet."
                                fi
                        fi
                break
                fi
        done
contigLength=`echo ${#newSeq}`
echo "The global add count is $(( $contigLength - $originalContigLength )), for a total of $contigLength nucleotides."
echo "This is contig $(($i + 1)) out of $scafListLength."
done
done
rm connected/empty; rm extended/empty; rm failed/empty; rm split_scaf/empty; rm temp_connected/empty
