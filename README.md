# 1D_KMC_matlab
A simple MATLAB script for solving a 1D time-dependent linearized phonon transport problem with the Kinetic Monte Carlo method

The system is a slab of Si material (dispersion and relaxation time data are located in the file dataSi.txt) 
of width L, with imposed tempeartures on both sides. The method makes use of the kinetic type algorithm
developed in :

<p><img src="../../Publications/whiteball.gif" alt="" align="bottom">
J-P. Péraud and N. G. Hadjiconstantinou,
<b>"An alternative approach to efficient simulation of micro/nanoscale phonon transport"</b>,
<em>Applied Physics Letters,</em> <b>101</b>, 153114, 2012.
</p>

MC1D.m is the main script. It uses the files select_mode.m and dataSi.txt, which must be in the MATLAB path 
To visualize a movie of the solution, uncomment the last 7 lines.

