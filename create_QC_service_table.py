import sqlite3
import pandas as pd

def list_tables(sql_cursor):
    sqlQ = "SELECT name FROM sqlite_master WHERE type='table'"
    sql_cursor.execute(sqlQ)
    tableList = [ t for l in cur.fetchall() for t in l ]
    return tableList

def create_service_table(sql_cursor):

    if 'service' in list_tables(sql_cursor):
        sql_cursor.execute("DROP TABLE service")
        
    sqlQ = (''' CREATE TABLE service 
                                (   procedure_id INTEGER PRIMARY KEY,
                                    date TEXT NOT NULL,
                                    type TEXT,
                                    is_pm TEXT,
                                    comment TEXT
                                )''')

    sql_cursor.execute(sqlQ)

    return True

def write_service(sql_cursor, df):
    
    sqlQ = ("INSERT INTO service VALUES (?, ?, ?, ?, ?)")
    
    for r in df.itertuples():
        print(list(r))
        r = [None,] + list(r)[1:]
        sql_cursor.execute(sqlQ, r)
    
    return True

if __name__ == '__main__':

    db_file = 'qc_qehf_st191107.db'

    conn = sqlite3.connect(db_file)
    conn.execute("PRAGMA foreign_keys = ON")
    cur = conn.cursor()
    
    create_service_table(cur)
    conn.commit()
    
##    df = pd.read_csv('Lumos_Cleaning.csv')
##    write_service(cur, df)
##    conn.commit()

    conn.close()

    
