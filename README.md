# qc-script-PD24
QC scripts for Proteome Discoverer 2.4

The scripts collect the information from the Proteome Discoverer search of a QC LC-MS run, save the key values into a database and create a visualization of the latest QC results. The database and visualizations are adapted for the Mascot searches and the injections of 50 ng of HeLa cell tryptic digest.

1) To start with, one needs to create an SQLite database. I have chosen to create one database file per actual mass spectrometer, but I can see why it can be handy to change the structure, put all instruments in the same file and table etc. The file should contain one or more tables for significantly different modes of use, for example the database for the Orbitrap Fusion Lumos mass spectrometer contains tables "lumos" and "lumos_faims", since these two modes are quite different and the QC results not very comparable between them. See the SQL columns below.

2) I have later added the "service" table to each file. The table was primarily created to hold infromation on the cleaning interventions, since they are important for the performance of a mass spectrometer.

3) One needs to specify the paths to the SQLite database files for each disctinct instrument mode, such as "fusion", "lumos", "lumos_faims" etc. For a local solution, one could simply place the database files on the network drive that is accessible throughout from any computer that needs access to it. The instrument name and mode of a QC run are inferred from the file name and form the information that is contained in the PD output.

4) The "reader_plotter" script is added to the Scripting Node, which is available in Proteome Discoverer since version 2.4:

<img src="https://github.com/dev-ev/QC_Script_PD2.4/blob/master/Screenshot_PD2.4_QC_ConsensusWF.PNG" alt="drawing" width="400"/>

Python 3.7 has been tested and used with the script. The script reads the output tables that have been saved as temporary files by PD, saves the key values to the SQLite database, and shows the main metrics from the current search as well as from the previous runs that correspond to the same instrument type. For convenience, the graphical report is saved as a png image, and the script launches the default image viewer program in order to show the report when it's ready:

<img src="https://github.com/dev-ev/QC_Script_PD2.4/blob/master/QC_graphical_report_example.png" alt="drawing" width="800"/>

Consider that only the key sums, averages and few other values are saved into the database, whereas the comprehensive information like the distribution of delta M vs retention time is only shown on the graphical report.

The main QC table contains the following columns:

    search_id INTEGER PRIMARY KEY,
    raw_file TEXT NOT NULL,
    file_date TEXT,
    search_date TEXT,
    instrument TEXT,
    protein_number INTEGER,
    peptide_number INTEGER NOT NULL,
    psm_number INTEGER NOT NULL,
    msms_number INTEGER NOT NULL,
    id_rate REAL,
    mean_psm_it_ms REAL,
    median_psm_it_ms REAL,
    mean_msms_it_ms REAL,
    median_msms_it_ms REAL,
    mean_mz_err_ppm REAL,
    median_mz_err_ppm REAL,
    mz_err_ppm_stdev REAL,
    total_prec_intensity REAL,
    mean_prec_intensity REAL,
    mean_sengine_score REAL,
    mean_peak_width REAL,
    peak_width_stdev REAL,
    pept_416 REAL,
    pept_425 REAL,
    pept_488 REAL,
    pept_495 REAL,
    pept_567 REAL,
    pept_652 REAL,
    pept_655 REAL,
    comment TEXT

The "service" table contains the following columns:

    procedure_id INTEGER PRIMARY KEY,
    date TEXT NOT NULL,
    type TEXT,
    is_pm TEXT,
    comment TEXT
