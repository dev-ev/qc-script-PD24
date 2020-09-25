# QC_Script_PD2.4
QC scripts for Proteome Discoverer 2.4

Created and tested by Egor Vorontsov. Questions and concerns by e-mail at gmail: yegor dot msu.

The scripts collect the information from the Proteome Discoverer search of a QC run, save the key values into a database and create a visualization of the latest QC results. The database and visualizations are adapted for the Mascot searches and the injections of 50 ng of Pierce HeLa digest

1) One needs to create an SQLite database first. We have chosen to create one database file per actual mass spectrometer. The file should contain one or more tables for significantly different modes of use, for example the database for the Orbitrap Fusion Lumos instrument contains tables "lumos" and "lumos_faims", since these two modes are quite differently and not very comparable between each other.

2) We have later added the "service" table to each file. The table was primarily created to hold infromation on the cleaning interventions, since they are important for the performance of the mass spectrometers.

3) One needs to specify the paths to the SQLite database files for each disctinct instrument mode, such as "fusion", "lumos", "lumos_faims" or others. We have simply placed the database files on the network drive that is accessible throughout the lab. The instrument name and mode of a QC run is inferred from the file name and form the information that is contained in the PD output.

4) The "reader_plotter" script is added to the Scripting Node, which is available in Proteome Discoverer since version 2.4. Python 3.7 has been working fine with the script. The script reads the output tables that has been saved as temporary files by PD, saves the key values to the SQLite database, and shows the main metrics from the current search as well as from the previous runs that correspond to the same instrument type in the databases. For convenience, the graphical report is saved as a png image, and the script launches the default image viewer program in order to show the report when it's ready. Consider that only the key sums, averages and other values are saved into the database, and the comprehensive information like the distribution of delta M vs retention time is only shown on the graphical report.
