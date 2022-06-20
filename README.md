# MACAROM

> MAximization of ClAsification of Raw Obfuscated Modulons

## Evaluation
The evaluation data set is taken from the RegulonDB for _Escherichia coli_ K12. This data set has been experimentally validated and can be found here: http://regulondb.ccg.unam.mx/external_data/MatrixAlignment/results/.

The BindingSite_Dataset.txt file contains the binding site sequence and position for each transcription factor.

The PSSM-Dataset-v4.0.txt file contains the PSSM matrix for each transcription factor.

## Reproduce MEME

* Download and install MEME (5.4.1).
* Fetch the sequences corresponding to an iModulon into a FASTA:

```bash
python scripts/to_fasta.py data/process/gene_seqs.csv INPUT_FASTA IMODULON

```

* Run

```bash
meme INPUT_FASTA -dna -oc OUTPUT_DIR -nostatus -time 14400 -mod anr -nmotifs 1 -minw 6 -maxw 50 -objfun classic -revcomp -markov_order 0
```
