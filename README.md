# Implementation of Smith-Waterman local alignment
This is a repository for a CompBio assignment for CSEP527
course assignment page: https://courses.cs.washington.edu/courses/csep527/20au/hw/hw2.html
The data was scraped from FASTA, and are located in `amino-acid-sequences` directory.


# How to use

- For string inputs, an example command usage is shown below
```
python src/smith_waterman.py --str-input -A "AKA" -B "ak" -o output/output.txt
```
Note: the input is not case sensitive.

- For file inputs, an example command usage is shown below
```
python src/smith_waterman.py --file-input -af amino-acid-sequences/P15172.fasta -bf amino-acid-sequences/Q10574.fasta -p -o output/output.txt
```
**Caution: `-p` flag calculates p-value from Fisher Yates shuffles, so for fasta file inputs or long string sequences, it may take up to 20 minutes**

- The shell script organizes the runs to generate HW2 report. In case you'd like to see more examples of the runs, please take a look at the shell script.