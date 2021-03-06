Genome Values Lookup

This algorithm screens a fasta-format file, looking for values in a given range, comparing it to the sum of values relative to the codons of a stretch of sequence with a given length.

Usage:
lookup.py [-h] -F <example.fasta> -O <output.csv> -C <CSC.csv> -L <10> --min <-0.5> --max <0.5>
-h, --help (optional): Show help message and exit
-F, --fasta (required): Fasta file to be analyzed
-O, --output (required): Output path in CSV format
-C, --csc (required): CSV-format file containing the CSC values for every codon
-L, --length (required): Length of the stretch to be screened for
--min (optional): The minimum value that will be counted as a match. If absent, the minimum will be -9999
--max (optional): The maximum value that will be counted as a match. If absent, the maximum will be 9999

example.fasta is an example file downloaded from https://www.yeastgenome.org/

Brown.csv, Cramer.csv, Gresham.csv and Young.csv are CSC example files obtained from the following mRNA half-life studies:
-Brown: Wang, Y., Wang, Y., Liu, C.L., Liu, C.L., Storey, J.D., Storey, J.D., Tibshirani, R.J., Tibshirani, R.J., Herschlag, D., Herschlag, D., Brown, P.O., Brown, P.O., 2002. Precision and functional specificity in mRNA decay. Proceedings of the National Academy of Sciences 99, 5860–5865. doi:10.1073/pnas.092538799
-Cramer: Sun, M., Sun, M., Schwalb, B., Schwalb, B., Schulz, D., Schulz, D., Pirkl, N., Pirkl, N., Etzold, S., Etzold, S., Larivière, L., Lariviere, L., Maier, K.C., Maier, K.C., Seizl, M., Seizl, M., Tresch, A., Tresch, A., Cramer, P., Cramer, P., 2012. Comparative dynamic transcriptome analysis (cDTA) reveals mutual feedback between mRNA synthesis and degradation. Genome Research 22, 1350–1359. doi:10.1101/gr.130161.111
-Gresham: Neymotin, B., Athanasiadou, R., Gresham, D., 2014. Determination of in vivo RNA kinetics using RATE-seq. RNA 20, 1645–1652. doi:10.1261/rna.045104.114
-Young: Holstege, F.C., Jennings, E.G., Wyrick, J.J., Lee, T.I., Hengartner, C.J., Green, M.R., Golub, T.R., Lander, E.S., Young, R.A., 1998. Dissecting the regulatory circuitry of a eukaryotic genome. Cell 95, 717–728.