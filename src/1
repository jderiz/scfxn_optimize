cpn=""
for item in "1" "4" "29(x3)" "15" "21(x4)"
do
        count=${item%%(*} # delete longest match of ( from end 18(x7) --> 18 
        #echo "count: ${count}"
        #echo "item: ${item}"
        if [[ $count != $item ]] # if there actually was something to delete
        then 
                #echo "no match"
                repeats=${item##*x} # delete longest match of x from start 18(x7) --> 7)
                repeats=${repeats%)*} # delete shortest match of ) from end 7) --> 7
                #echo "repeats: ${repeats}"
                for i in $(eval echo "{1..$repeats}") # make iterator list
                do
                #        echo "cpn: ${cpn}"
                        cpn="${cpn} ${count}" # append count to result string
                done
        else
        cpn="${cpn} ${count}"

        fi
done
echo "${cpn}"
