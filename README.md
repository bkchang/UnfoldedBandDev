# UnfoldedBandDev
This script computes the average root-mean square deviation per eigenvalue of a BandUP output file from a primitive cell (PC) spectrum.

1. Choose an k-mesh for your PC calculation and calculate it with VASP.
2. Use BandUP to unfold a supercell (SC) band structure onto the k-points of the k-mesh.
3. The unfolded spectral data are in the output file "BandUP/step-4/plot/../...symmetry_averaged....dat".
4. Copy the IBZKPT and EIGENVAL of PC calculation, the unfolded SC spectral data, and this python script in the same directory.
5. Check that the parameters in the head of the python script are set to your desired setting.
6. Run the script.

# Note
* This script omits degeneracies in PC spectrum. If one needs it, one should implement it by him/herself.
* The output figures are as follows:
...1. PC.png --- PC spectrum
...2. SC_raw.png --- SC intensity profile before preprocessing
...3. SC_refined.png --- Preprocessed SC intensity profile with noise removal and bottom of valence band alignment
