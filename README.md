DISCLAIMER: this program was created as a school assignment at FEE CTU Prague, within the Bioinformatics course

AUTHOR:     Jiří Richter

  This Ruby script enables You to globally and locally align 2 sequences of DNA.
  Global alignment uses the Needleman-Wunsch algorithm.
  Local alignment uses the Smith-Waterman algorithm.
  
  Input: 
  
    - match, mismatch and gap coefficients
    - two DNA sequences (GATCgatc only)
    
  Output:  
  
    - matrix with filled scoring values
    - path array
    - alignment of the two sequences as arrays


Sample output
```
Input the MATCH scoring coefficient>
3
Input the MISMATCH scoring coefficient>
-2
Input the GAP scoring coefficient>
-5
Input first DNA sequence>
gattaca
Input second DNA sequence>
ata
The score matrix:
    -    A    T    A
-   0   -5  -10  -15 
A  -5    3   -4   -6 
C -10   -4    1   -6 
A -15   -6   -6    4 
T -20  -13   -3   -3 
T -25  -20   -5   -5 
A -30  -22  -12   -2 
G -35  -29  -19   -9 


Path (with upward preference):
[-9, -2, -5, -3, -6, -4, 3, 0]

DNA alignment:
["-", "A", "-", "T", "-", "-", "A"]
["G", "A", "T", "T", "A", "C", "A"]

Score: -26
The score matrix:
    -    A    T    A
-   0    0    0    0 
A   0    3    0    0 
C   0    0    1    0 
A   0    0    0    4 
T   0    0    0    0 
T   0    0    0    0 
A   0    0    0    0 
G   0    0    0    0 

Path (with upward preference):
[1, 3, 0]

DNA alignment:
["A", "T", "A"]
["A", "C", "A"]
```
