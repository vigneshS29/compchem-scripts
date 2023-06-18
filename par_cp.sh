#!/bin/bash  

function traverse () {
for file in "$1"/*; do
        if [ -d "$file" ]; then
                if [ ! -d "$2"/"$(echo $file | rev | cut -d '/' -f1 | rev)" ]; then
                        mkdir  "$2"/"$(echo $file | rev | cut -d '/' -f1 | rev)"
                fi
                traverse "$file" "$2"/"$(echo $file | rev | cut -d '/' -f1 | rev)" "$3"
        else                                                                                                                                                                                                                     
                cp "$file" "$2" &
		pids+=($!)
		counter=$((counter+1))
                counter_tot=$((counter_tot+1))
                
                if [ "$counter" == "$3" ]; then
                        for pid in "${pids[@]}"; do
                                wait "$pid"
                        done
                        pids=()
                        counter=0 
                fi          
        fi
done
}

function main() {
pids=()
counter=0
counter_tot=0
if [ ! -d "$2" ]; then
        mkdir ./"$2"
fi
traverse "$1" "$2" "$3"
echo "Copied" $counter_tot "files"
}

main "$1" "$2" "$3"
