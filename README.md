# Hamiltonian truncation tutorial for the CERN workshop (May 2022)

This is an example Hamiltonian truncation script. It applies to both the 2d massive scalar field on R x [0,L] and its counterpart on AdS_2.
The script computes two spectral densities that appear at second order in perturbation theory.
To run this code, you need to have Python 3 installed along with the packages NumPy and Matplotlib.

To run the code, clone the repository (or download the files to a directory) and type

    python3 run.py
  
in a terminal.

Several parameters can be changed: the bare mass \mu, the cutoff \Lambda and the geometry.
These settings can be changed in the top section of run.py.
To change the geometry, set either

    system = "strip"
  
or 

    system = "AdS"
  
in the top of the file run.py. 

Finally, you will notice that the AdS mass gap computed using the spectral density gives a wrong result.
In our paper, we explain how to fix this. To implement the prescription from the paper, locate the call

    (...) = buildSpectralDensities(mu,cutoff,system,shift = False)

inside run.py, and change it to

    (...) = buildSpectralDensities(mu,cutoff,system,shift = True)

You should be able to see the difference at the level of the plot of the spectral densities.
