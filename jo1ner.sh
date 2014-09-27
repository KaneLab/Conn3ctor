#This script will iterate through a "connected" directory made by connector.
#Joining together scaffolds that match query.
#Salt and Pepper included (error checking).

#Input is a file from connected folder

if [ `grep -c joined $1` -ne 0 ];then
        echo "Infile $1 already joined.  Moving on."
        break
fi
#number of headers
headCount=`grep -c ">" $1`

#get header name
header=`head -n1 $1 | tr -d ">"`

#get the query sequences that is the match between two contigs
seq=`grep "The query that matches both contigs is" $1 | awk -F ": " '{print $2}'`
#ensure the query is only found three times
matches=`grep -o $seq $1 | wc -l`

#initialize product length
prodLen=0

#print error if query is found more than three times
if [ $matches -ne 3 ] || [ $headCount -ne 2 ];then
        echo "Infile $1: Either the query was found more than three times, or the number of contigs is $headCount.  Fuck this, I'm out."
else

        #get length of the two contigs
        c1=`grep -m 1 -A1 ">" $1 | tail -n1`
        c1Length=${#c1}

        c2=`grep -A1 ">" $1 | tail -n1`
        c2Length=${#c2}

        #get max of the two contigs
        max=0
        if [ $c1Length -gt $c2Length ];then
                max=$c1Length
                maxHead=`grep -m 1 -A1 ">" $1 | head -n1 | tr -d ">"`
                grep -m 1 -A1 ">" $1 > ${maxHead}_contig.txt
        else
                max=$c2Length
                maxHead=`grep -A1 ">" $1 | tail -n2 | head -n1 | tr -d ">"`
                grep -A1 ">" $1 > ${maxHead}_contig.txt
        fi

        sed "s/$seq/\n$seq\n/g" $1 > ${header}_joined_temp.txt
        head -n2 ${header}_joined_temp.txt > ${header}_joined_top.txt
        tail -n2 ${header}_joined_temp.txt > ${header}_joined_bottom.txt
        cat ${header}_joined_top.txt ${header}_joined_bottom.txt  | tr -d "\n\t^M " | sed "s/$header/$header\n/g" > ${header}_joined.txt
        rm  ${header}_joined_temp.txt; rm  ${header}_joined_top.txt; rm ${header}_joined_bottom.txt

        product=`tail -n1 ${header}_joined.txt`
        prodLen=${#product}

        if [ $prodLen -lt $max ];then
                echo "Infile $1: The join is smaller than the largest contig! Moving the joined file into the badConnections directory.  A new file called ${maxHead}_contig.txt has been made because this is your largest contig"
                mv ${maxHead}_contig.txt > ./extended/${maxHead}_contig.txt
                #check for existence of badConnections directory
                check=`ls -d badConnections | wc -l`
                if [ $check -ne 0 ];then
                        mv $1 badConnections/
                else
                        mkdir badConnections/
                        mv $1 badConnections/
                fi

        else
                rm ${maxHead}_contig.txt
                echo "Infile $1:  Success!  Moving $1 into infile directory for your records."
                iCheck=`ls -d infile | wc -l`
                if [ $iCheck -ne 0 ];then
                        mv $1 infile/
                else
                        mkdir infile/
                        mv $1 infile/
                fi


        fi
fi
