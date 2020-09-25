import csv
import json
import os
import pandas as pd
import warnings

class InputTables:
    ''' Reads the info about the input tables.
        Reads the tables as Pandas dataframes
        Stores and returns the tables as a dictionary, or one by one
    '''

    def __init__(self, table_list):
        self.df_dict = {}
        self.table_list_of_dicts = table_list
        for table_description in self.table_list_of_dicts:
            table_name = table_description['TableName']
            df = pd.read_csv(table_description['DataFile'],delimiter='\t')
            self.df_dict[table_name] = df

    def print_all_dfs(self):
        for i in self.df_dict.keys():
            print('Table Name', i)
            print('Size', self.df_dict[i].shape)
##            print(self.df_dict[i].columns)
            print(self.df_dict[i])
            print('Columns:', self.df_dict[i].columns)
##            for cn in self.df_dict[i].columns:
##                print(self.df_dict[i][cn])

    def return_all_tables(self):
        return self.df_dict

    def save_all_tables(self):
        for table_description in self.table_list_of_dicts:
            tn = table_description['TableName']
            basename = table_description['DataFile']
            basename = basename.split('/')[-1][:-4]
            fn = basename + '_test_out.xlsx'
            self.df_dict[tn].to_excel(fn)
            print('Tables saved successfully as Excel files!')
        return True

class NodeArgs:
    ''' Reads json input file from the PD
        Stores and returns the info from the json input tables
    '''

    def __init__(self, json_fname):
        self.loading, self.jdict = self.load_json_file(json_fname)
        if self.loading is False:
            warnings.warn('Json dictionary is empty',
                          UserWarning)
        if self.loading is True:
            self.parse_in_json_dict()

    def load_json_file(self, json_fname):
        jobj =  ''
        try:
            jfile = open(json_fname,'r')
            try:
                jobj = json.load(jfile)
            except:
                warnings.warn('Could not read the object as json',
                              UserWarning)
        except:
            warnings.warn(('Could not read the file '+json_fname),
                          UserWarning)
        if jobj != '':
            return (True,jobj)
        else:
            return (False,None)

    def loading_successful(self):
        return self.loading

    def parse_in_json_dict(self):
        self.cur_wf_id = self.jdict['CurrentWorkflowID']
        self.responce_path = self.jdict['ExpectedResponsePath']
        self.result_path = self.jdict['ResultFilePath']
        self.node_params = self.jdict['NodeParameters']
        self.ver_from_injson = self.jdict['Version']
        self.intable_list = self.jdict['Tables']
        return True

    def return_all_table_properties(self):
        return self.intable_list

    def save_json_input_dict(self):
        efn = 'E:\\Projects\\Thermo_Script_Integration\\PD24\\Custom_Scipts\\json_input_dictionary.txt'
        with open(efn, 'w') as fh:
            fh.write('Json input dict was:\n')
            for akey in self.jdict.keys():
                fh.write(akey + ':\n')
                fh.write(str(self.jdict[akey]) + '\n')
            fh.write('Separate table listing:\n')
            for i in self.intable_list:
                fh.write('Table <'+ str(i['TableName']) + '>:\n')
                fh.write(str(i) + '\n')
        print('Json input dictionary saved as', efn)
        return True

def testing_load_example_files():
    df_dict = {}
    df_dict['Proteins'] = pd.read_csv('TargetProtein.txt',delimiter='\t')
    df_dict['Peptide Groups'] = pd.read_csv('TargetPeptideGroup.txt',delimiter='\t')
    df_dict['PSMs'] = pd.read_csv('TargetPeptideSpectrumMatch.txt',delimiter='\t')
    df_dict['MS/MS Spectrum Info'] = pd.read_csv('MSnSpectrumInfo.txt',delimiter='\t')
    df_dict['Input Files'] = pd.read_csv('WorkflowInputFile.txt',delimiter='\t')
    
    return df_dict

if __name__ == '__main__':

    a = NodeArgs('node_args_psm.json')
    if a.loading_successful() is True:
        b = InputTables(a.return_all_table_properties())
        b.print_all_dfs()
##    dfs = testing_load_example_files()
##    for i in dfs.keys():
##        print(i)
##        print(dfs[i].columns)
##        print(dfs[i].shape)
