# shell script for the runs

rm -rf output/output.txt && touch output/output.txt

# III.a
echo "III.a" >> output/output.txt
python src/smith_waterman.py --str-input -B "ddgearlyk" -A "deadly" -p -o output/output.txt

# III.b
echo "\n\n\nIII.b" >> output/output.txt
files=( P15172 P17542 P10085 P16075 P13904 Q90477 Q8IU24 P22816 Q10574 O95363 )

for i in $files; do
  for j in $files; do
    if [[ $i != $j ]]; then
      python src/smith_waterman.py --file-input -af amino-acid-sequences/$i.fasta -bf amino-acid-sequences/$j.fasta -o output/output.txt
    fi
  done
done

# III.c
echo "\n\n\nIII.c" >> output/output.txt
python src/smith_waterman.py --file-input -af amino-acid-sequences/P15172.fasta -bf amino-acid-sequences/Q10574.fasta -p
python src/smith_waterman.py --file-input -af amino-acid-sequences/P15172.fasta -bf amino-acid-sequences/O95363.fasta -p