# icecores_diffusion
Diffusion-correction code for water isotopes in ice cores, as well as extrema picking for seasonally-resolved data.

Diffusion-correction code developed by Sigfus Johnsen, 
University of Copenhagen. which uses Maximum Entropy 
Methods (MEM) to invert an observed power density spectrum.
Code provided in 2011 by Bo M. Vinther, University of Copenhagen.
Code adapted in 2012 by Tyler R. Jones, University of Colorado,
to include extrema picking.

Used in the following manuscript to reconstruct seasonal
isotopic-variability (maximum summer and minimum winter values)
in the WAIS Divide ice core.

Jones, T. R., Cuffey, K. M., Roberts, W. H. G., Markle, B. R., 
Steig, E. J., Stevens, Cd18O_original.txt. M., Valdes, P. J., Fudge, T. J., 
Sigl, M., Hughes, A. G., Morris, V., Vaughn, B. H., Garland, J., 
Vinther, B. M., Rozmiarek, K. S., Brashear, C. A., & 
White, J. W. C. (2022, in press). Seasonal temperatures in 
West Antarctica during the Holocene. Nature.

Jones, T. R. (2022) "Seasonal temperatures in West Antarctica 
during the Holocene " U.S. Antarctic Program (USAP) Data Center. 
doi: https://doi.org/10.15784/601603.

Input files include d18O_diffused.txt, to be diffusion-corrected,
and d18O_original.txt, a hypothetical pre-diffusion signal for
comparison purposes to the diffusion-correction and as an 
example for understanding fitting variables (length of 
deconvolution filter, cut-off frequency, and diffusion length)
within a window of data.

In ice core science, cut-off frequency is determined by 
quantifying the frequency at which instrumental noise 
overwhelms the climate signal, and then choosing a sufficiently 
lower frequency to avoid this noise. The diffusion length
is quantified by the methods in Jones et al. 2017b (see below).
The length of deconvolution filter is usually set to 100.

Jones, T. R., Cuffey, K. M., White, J. W. C., Steig, E. J., 
Buizert, C., Markle, B. R., McConnell, J. R. & Sigl, M. (2017b). 
Water isotope diffusion in the WAIS Divide ice core during 
the Holocene and last glacial, J. Geophys. Res. Earth Surf., 122.
