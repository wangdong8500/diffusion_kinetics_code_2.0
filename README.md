# diffusion_kinetics_code_2.0
This is a code for sovling one dimensional diffusion kinetics equations with application to Jovian planets atmosphere chemistry.
To test this module, go to the 'work' directory and run $python main.py. 
Two input files are needed. The PT_Jupiter.csv is a one dimensional pressure temperature grid of Jupiter's atmosphere. The mechanism.cti is the chemistry input file needed by cantera. Cantera must be installed to properly run this module. The output files are written into 'output' directory. 
Other chemistry can be accounted by changing the *.cti file. Other atmospheres can be accounted by changing the *.csv file. More flexible output can be done by editing the 'main.py' file

