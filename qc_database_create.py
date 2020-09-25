import sqlite3

def create_tables(sql_cursor, table_names):

    for table in table_names:
        table_command = (''' CREATE TABLE IF NOT EXISTS ''' + table +
                         ''' (
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
                                    )''')

        sql_cursor.execute(table_command)
        populate_with_simulated_vals(sql_cursor, table)
    return True

def populate_with_simulated_vals(sql_cursor, table_name):

    sql_command = ("INSERT INTO " + table_name +
                   " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)")
    
    d = '2019/10/24'
    s = '2019/11/07'
    i = 'Virtual HF'
    simulated_vals_tuple = ((None,'Nothing_01.raw',d,s,i,2000,9000,9500,25000,9500/25000
                             ,20,18,25,22,-3.2,-2.2,1,2500000,65685,25.5,20,5.5,
                             28,24,19,35,16,32,40,'Simulated'),
                            (None,'Nothing_02.raw',d,s,i,2000,10000,10000,30000,10000/30000
                             ,18.2,17.5,20.1,19.5,-1.2,-0.8,0.26,2800000,80000,26.2,21.1,6.8,
                             28,24,19,35,16,32,40,'Simulated'),
                            (None,'Nothing_03.raw',d,s,i,1500,6000,6200,22000,6200/22000
                             ,35.1,32.6,43.1,35,1.6,1.5,0.6,2900000,91000,28.1,20.3,7.1,
                             28,24,19,35,16,32,40,'Simulated'),
                            (None,'Nothing_04.raw',d,s,i,1826,7586,8123,29365,8123/29365
                             ,20,19,23,18.6,0.78,0.65,0.1,30000000,100020,25.17,18.6,5.5,
                             28,24,19,35,16,32,40,'Simulated'),
                            (None,'Nothing_05.raw',d,s,i,2230,11258,12000,32000,12000/32000
                             ,13,12.8,15.6,15.6,4.6,4.2,1.6,3100000,120055,28.42,19.35,4.78,
                             28,24,19,35,16,32,40,'Simulated'))

    for val_tuple in simulated_vals_tuple:
        sql_cursor.execute(sql_command, val_tuple)
        
    return True

if __name__ == '__main__':

##    db_file = 'qc_lumos_st191029.db'
##    tables = ('lumos','lumos_faims')
    db_file = 'qc_qehf_st191107.db'
    tables = ('qehf',)

    conn = sqlite3.connect(db_file)
    cur = conn.cursor()

    create_tables(cur, tables)
    
    conn.commit()
    conn.close()

    
