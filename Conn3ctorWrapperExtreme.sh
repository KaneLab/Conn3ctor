#Inputs are:

# $1=parameter file
# $2=genome file with blocks of N removed
# $3=number of iterations to run Conn3ctor

#change input filename to work with wrapper loop
cp $2 ${2}_0.tempGenome

#start wrapper loop
for (( i=0; i<$3; i++));do
        bash Conn3ctor.sh $1 ${2}_${i}.tempGenome
        echo "y"

        #implement jo1ner on each file in ./connected/
        for j in `ls connected`;do
                bash jo1ner.sh connected/$j
                mv ${j}_joined.txt connected/
        done

        #get next value of i
        iPlus=$(( $i + 1 ))

        #gather all of the
        cat connected/*joined* > ${2}_${iPlus}.tempGenome
        cat extended/* >> ${2}_${iPlus}.tempGenome
        cat failed/* >> ${2}_${iPlus}.tempGenome

        #insert carriage returns between contigs
        sed -i 's/>/\n>/g' ${2}_${iPlus}.tempGenome
done
