medianNums=`grep -v '>' $1 | sed 's/\s/\n/g' | sort -g` #get rid of header line in the cap3 quality file and put each number on a separate line for sorting
medianArray=($medianNums) #cast into an array
medLength=${#medianArray[@]} #get array length
medpos=$((medLength/2)) #determine position of median value
median=${medianArray[$medpos]} #extract median value
echo $median
