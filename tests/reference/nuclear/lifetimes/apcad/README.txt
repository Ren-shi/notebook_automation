
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~ Example application for APCAD ~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two example fits to simulated spectra are provided in the files

	OneShape_example.dsam and
	ThreeShapes_example.dsam

for demonstration of the program APCAD. These files can be
opened via the "Load Session" button in APCAD and provide a
completed fit to simulated spectra.


The simulated spectra are located in the subfolder 'Input'.
They were calculated with the following parameters:

Spectrum 'OneShape_example.spc'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Transition energy: 	1000keV
Level lifetime:		1ps
Number of events: 	7400
Ang. distr. coeff. A20: 0.4
Ang. distr. coeff. A20: -0.8
Detector response FWHM: 3keV
Background function:    400 - 0.252*x + 5.6e-5*x²

Spectrum 'ThreeShapes_example.spc'
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
as before, but with 2 contaminants:

--First contaminant--
Transition energy: 	1025keV
Level lifetime:		0.5ps
Number of events: 	5000
Ang. distr. coeff. A20: 0.3
Ang. distr. coeff. A20: 0.0

--Second contaminant--
Transition energy: 	985keV
Level lifetime:		1.6ps
Number of events: 	6000
Ang. distr. coeff. A20: 0.2
Ang. distr. coeff. A20: -0.15


Stopping matrix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The stopping matrix in contained in the file 'Stopmat_example.root'
in the subfolder "Stop".

To calculate the Stopping Matrix by yourself, create a symbolic
link to the StopSim program in the subfolder "Stop" and run StopSim
by the command

	./StopSim Macro_StopSim_example.mac

This re-creates the file 'Stopmat_example.root'. The excitation cross
section is read from the file 'XSec_example.dat'. 


Projection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The projected stopping matrix is contained in the 
file 'Projection_example.proj' in the subfolder "Stop".

To perform the projection by yourself, you have make sure the 
detector geometry is set. For the example data, a position-sensitive
detector covering an azimuthal angular from -15° to 15° in each 
angular bin. A detector definition file 'Detector_example.det' is 
located in the subfolder "Input".

