# Implementation of Smith-Waterman local alignment
This is a repository for a CompBio assignment for CSEP527
course assignment page: https://courses.cs.washington.edu/courses/csep527/20au/hw/hw2.html

# How to use

For string inputs, an example command usage is shown below
```
python src/smith_waterman.py --str-input -A "AKA" -B "ak"
```

For file inputs, an example command usage is shown below
```
python src/smith_waterman.py --file-input -af amino-acid-sequences/P15172.fasta -bf amino-acid-sequences/Q10574.fasta -p
```
**Caution: `-p` flag calculates p-value from Fisher Yates shuffles, so for fasta file inputs or long string sequences, it may take up to 20 minutes**