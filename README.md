# Fasta Cleaner

Cleans fasta files.
All sequences of newlines and carriage returns are replaced by single newlines.
The Record headers are left unchanged, and the sequences are transformed into upper case, and all characters that are not A, C, G or T are removed.

While characters are removed, the line width of the input file is left intact.
It is guessed from the width of the first input sequence line, and all subsequent sequence strings are adjusted accordingly.
The adjustment happens via moving line breaks, and not via removing valid sequence characters.

## Example

Input:

```plain
\r>WGCaC\n\nAACCcxXAA\naacc\n.ef34\nCGG\ntgtcgcgtagcgtgatcgtgtagtcgtag\r.\r>f\nTTT
```

Output

````plain
>WGCaC\nAACCCAAAA\nCCCGGTGTC\nGCGTAGCGT\nGATCGTGTA\nGTCGTAG\n>f\nTTT\n
```

## Known Issues

If the first sequence line is shorter than the line width of the input fasta file, then the sequence lines in the output fasta file will be adjusted accordingly.
