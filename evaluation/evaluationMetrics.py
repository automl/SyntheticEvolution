import os
import sys
import ast

ROOT = os.path.dirname(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from plots.plotCreator import PlotCreator
from database.databaseMethods import DatabaseMethods
import matplotlib.pyplot as plt
import pandas as pd

def max_length(data):
    data_max_length = []
    for exp_db_id, x, rna_length in data:
        if isinstance(rna_length, str) and rna_length.startswith('[') and rna_length.endswith(']'):
            rna_length = max(ast.literal_eval(rna_length))
        data_max_length.append((exp_db_id, x, rna_length))
    return data_max_length

class ProteinMetrics(DatabaseMethods, PlotCreator):
    def __init__(self, id: str = "", all_rmsd: str = "", protein_rmsd: str = "", protein_lddt: float = None,
                 protein_tm_score: float = None,  **kwargs):
        DatabaseMethods.__init__(self)
        if 'msa_option' in kwargs or 'single_chain_only' in kwargs:
            PlotCreator.__init__(self, 'evaluation_metrics',
                                 kwargs.get('msa_option'),
                                 kwargs.get('single_chain_only'))
        self.id = id
        self.all_rmsd = all_rmsd
        self.protein_rmsd = protein_rmsd
        self.protein_lddt = protein_lddt
        self.protein_tm_score = protein_tm_score

    def load_from_db(self, select_table):
        columns = ["Complex_RMSD", "Protein_RMSD", "Protein_LDDT", "Protein_TM"]
        condition = "exp_db_id = ?"
        params = (self.id,)
        data = self.get_table_columns(select_table, columns, condition, params)

        if data:
            (self.all_rmsd, self.protein_rmsd, self.protein_lddt, self.protein_tm_score) = data

    def get_protein_rmsd_values(self, select_table):
        # rows = self.get_table_columns(select_table, ["Protein_RMSD"], condition=None)
        rows = self.get_filtered_values(select_table, "Protein_RMSD")
        protein_rmsd = [
            value_list[0] for row in rows if row[0]
            for value_list in [ast.literal_eval(row[0])]
            if isinstance(value_list, list) and value_list
        ]
        return protein_rmsd

    def get_protein_lddt_values(self, select_table):
        # rows = self.get_table_columns(select_table, ["Protein_LDDT"], condition=None)
        rows = self.get_filtered_values(select_table, "Protein_LDDT")
        lddt_values = [row[0] for row in rows if row[0] not in (None, '')]
        return lddt_values

    def get_protein_tm_values(self, select_table):
        # rows = self.get_table_columns(select_table, ["Protein_TM"], condition=None)
        rows = self.get_filtered_values(select_table, "Protein_TM")
        protein_tm_score = [row[0] for row in rows if row[0] not in (None, '')]
        return protein_tm_score

class RNAMetrics(DatabaseMethods, PlotCreator):
    def __init__(self, id: str = "", all_rmsd: str = "", rna_rmsd: str = "", rna_lddt: float = None,
                 rna_tm_score: float = None, rna_inf: str = "", rna_di: float = None, **kwargs):
        DatabaseMethods.__init__(self)
        if 'msa_option' in kwargs or 'single_chain_only' in kwargs:
            PlotCreator.__init__(self, 'evaluation_metrics',
                                 kwargs.get('msa_option'),
                                 kwargs.get('single_chain_only'))
        self.id = id
        self.all_rmsd = all_rmsd
        self.rna_rmsd = rna_rmsd
        self.rna_lddt = rna_lddt
        self.rna_tm_score = rna_tm_score
        self.rna_inf = rna_inf
        self.rna_di = rna_di

    def load_from_db(self, select_table):
        columns = ["Complex_RMSD", "RNA_RMSD", "RNA_LDDT", "RNA_TM", "RNA_INF", "RNA_DI"]
        condition = "exp_db_id = ?"
        params = (self.id,)
        data = self.get_table_columns(select_table, columns, condition, params)

        if data:
            (
                self.all_rmsd,
                self.rna_rmsd,
                self.rna_lddt,
                self.rna_tm_score,
                self.rna_inf,
                self.rna_di
            ) = data

    def get_rna_rmsd_values(self, select_table):
        rows = self.get_filtered_values(select_table, "RNA_RMSD")
        rna_rmsd = [
            value_list[0] for row in rows if row[0]  # row is a tuple (value,)
            for value_list in [ast.literal_eval(row[0])]
            if isinstance(value_list, list) and value_list
        ]
        return rna_rmsd

    def get_rna_rmsd(self):
        return self.rna_rmsd

    def get_rna_lddt_values(self, select_table):
        # rows = self.get_table_columns(select_table, ["RNA_LDDT"], condition=None)
        rows = self.get_filtered_values(select_table, "RNA_LDDT")
        lddt_values = [row[0] for row in rows if row[0] not in (None, '')]
        return lddt_values

    def get_rna_tm_values(self, select_table):
        # rows = self.get_table_columns(select_table, ["RNA_TM"], condition=None)
        rows = self.get_filtered_values(select_table, "RNA_TM")
        rna_tm_score = [row[0] for row in rows if row[0] not in (None, '')]
        return rna_tm_score

    def get_x_and_y(self, x, y, select_table):
        # print("@£$%@£%@", x, y)
        condition = f"{x} IS NOT NULL AND {y} IS NOT NULL"
        data = self.get_filtered_values(
            select_table,
            ["exp_db_id", y, x],
            base_condition=condition
        )
        return max_length(data)

    def get_x_and_y_condition(self, x, y, condition, select_table):
        # print("@£$%@£%@", x, y)
        condition = f"{condition} AND {x} IS NOT NULL AND {y} IS NOT NULL"
        data = self.get_filtered_values(
            select_table,
            ["exp_db_id", y, x],
            base_condition=condition
        )
        return max_length(data)

    def get_rna_inf_values(self, select_table):
        rows = self.get_filtered_values(select_table, "RNA_INF")
        # rows = self.get_table_columns(select_table, ["RNA_INF"], condition=None)
        all_inf = []
        wc_inf = []
        nonWC_inf = []
        stacking_inf = []
        for row in rows:
            if row[0]:
                try:
                    value_list = ast.literal_eval(row[0])
                    if isinstance(value_list, list) and len(value_list) > 0:
                        all_inf.append(value_list[0])
                        wc_inf.append(value_list[1])
                        nonWC_inf.append(value_list[2])
                        stacking_inf.append(value_list[3])
                except (ValueError, SyntaxError):
                    continue
            else:
                continue
        return all_inf, wc_inf, nonWC_inf, stacking_inf

    def get_rna_inf_values_with_extra(self, select_table, extra_column):
        # query = f"SELECT exp_db_id, {extra_column}, RNA_INF FROM unified_pred WHERE SourceTable = ?"
        # rows = self.get_table_columns(select_table, ["exp_db_id", extra_column, "RNA_INF"], condition=None)
        rows = self.get_filtered_values(select_table, ["exp_db_id", extra_column, "RNA_INF"])

        exp_db_ids = []
        x_values = []  # This will be 'extra_value'
        y_values_all_inf = []  # This will be 'RNA_INF'
        wc_inf = []
        nonWC_inf = []
        stacking_inf = []

        for row in rows:
            exp_db_id = row[0]
            extra_value = row[1]
            rna_inf_value = row[2]

            # Check if both extra_value and rna_inf_value are present
            if extra_value is not None and rna_inf_value:
                try:
                    value_list = ast.literal_eval(rna_inf_value)
                    if isinstance(value_list, list) and len(value_list) > 0:
                        # Append to lists only if valid values are found
                        exp_db_ids.append(exp_db_id)
                        x_values.append(extra_value)
                        y_values_all_inf.append(value_list[0])  # y-values for plotting
                        wc_inf.append(value_list[1])
                        nonWC_inf.append(value_list[2])
                        stacking_inf.append(value_list[3])
                except (ValueError, SyntaxError):
                    continue

        # Ensure all lists are the same length
        min_length = min(len(x_values), len(y_values_all_inf), len(wc_inf), len(nonWC_inf), len(stacking_inf))
        return (exp_db_ids[:min_length], x_values[:min_length], y_values_all_inf[:min_length],
                wc_inf[:min_length], nonWC_inf[:min_length], stacking_inf[:min_length])

    def get_rna_di_values(self, select_table):
        # rows = self.get_table_columns(select_table, ["RNA_DI"], condition=None)
        rows = self.get_filtered_values(select_table, "RNA_DI")
        rna_di = [row[0] for row in rows if row[0] not in (None, '')]
        return rna_di

    def get_di_and_length(self, select_table):
        condition = "RNA_DI IS NOT NULL AND RNA_DI != '' AND RNALength IS NOT NULL"
        data = self.get_filtered_values(
            select_table, 
            ["exp_db_id", "RNA_DI", "RNALength"], 
            base_condition=condition
        )
        return max_length(data)

    def get_x_and_length(self, x, select_table):
        condition = f"{x} IS NOT NULL AND {x} != '' AND RNALength IS NOT NULL"
        data = self.get_filtered_values(
            select_table, 
            ["exp_db_id", x, "RNALength"], 
            base_condition=condition
        )
        return max_length(data)

    def get_rmsd_and_length(self, select_table):
        condition = "RNA_RMSD IS NOT NULL AND RNA_RMSD != '' AND RNALength IS NOT NULL"
        rows = self.get_filtered_values(
            select_table, 
            ["exp_db_id", "RNA_RMSD", "RNALength"], 
            base_condition=condition
        )
        processed_data = []
        for row in rows:
            exp_db_id = row[0]
            rna_rmsd_str = row[1]
            rna_length = row[2]
            try:
                rna_rmsd = float(rna_rmsd_str.strip('[]').split(',')[0]) #pymol: 0, tmscore:1
                processed_data.append((exp_db_id, rna_rmsd, rna_length))
            except (ValueError, IndexError):
                continue
        return max_length(processed_data)

    def get_msa_rmsd_comparison(self, table):
        """Get RMSD values grouped by MSA status"""
        query = f"""
        SELECT RNA_RMSD, af3_rna_MSA 
        FROM {table}
        WHERE RNA_RMSD IS NOT NULL AND af3_rna_MSA IS NOT NULL
        """
        rows = self.execute_query(query)
        
        msa_rmsd = {'+MSA': [], '-MSA': []}
        for rmsd, msa in rows:
            if rmsd and msa:
                try:
                    # Convert string representations to lists if necessary
                    if isinstance(rmsd, str):
                        rmsd_values = eval(rmsd)
                        rmsd = max(rmsd_values) if isinstance(rmsd_values, list) else float(rmsd)
                    if isinstance(msa, str):
                        msa_values = eval(msa)
                        # Classify as -MSA if any value is 0, otherwise +MSA
                        is_msa = 0 if (isinstance(msa_values, list) and 0 in msa_values) else 1
                    else:
                        is_msa = int(msa)
                    
                    # Add to appropriate list
                    if is_msa == 1:
                        msa_rmsd['+MSA'].append(float(rmsd))
                    else:
                        msa_rmsd['-MSA'].append(float(rmsd))
                except (ValueError, SyntaxError) as e:
                    print(f"Error processing values - RMSD: {rmsd}, MSA: {msa}")
                    continue
        
        return msa_rmsd['+MSA'], msa_rmsd['-MSA']

    def get_msa_metric_comparison(self, table, metric_column):
        """Get metric values grouped by MSA status"""
        query = f"""
        SELECT {metric_column}, af3_rna_MSA 
        FROM {table}
        WHERE {metric_column} IS NOT NULL AND af3_rna_MSA IS NOT NULL
        """
        rows = self.execute_query(query)
        
        metric_values = {'+MSA': [], '-MSA': []}
        for value, msa in rows:
            if value and msa:
                try:
                    # Handle RMSD values differently (they're stored as lists)
                    if metric_column == 'RNA_RMSD':
                        if isinstance(value, str):
                            value_list = eval(value)
                            if isinstance(value_list, list) and value_list:
                                value = float(value_list[0])  # Take first value for RMSD
                            else:
                                value = float(value)
                        else:
                            # For other metrics (TM, LDDT)
                            if isinstance(value, str):
                                value = float(value)
                    
                    # Process MSA values
                    if isinstance(msa, str):
                        msa_values = eval(msa)
                        # Classify as -MSA if any value is 0, otherwise +MSA
                        is_msa = 0 if (isinstance(msa_values, list) and 0 in msa_values) else 1
                    else:
                        is_msa = int(msa)
                    
                    # Add to appropriate list
                    if is_msa == 1:
                        metric_values['+MSA'].append(float(value))
                    else:
                        metric_values['-MSA'].append(float(value))
                except (ValueError, SyntaxError) as e:
                    print(f"Error processing values - {metric_column}: {value}, MSA: {msa}")
                    continue
        
        return metric_values['+MSA'], metric_values['-MSA']

    def get_msa_inf_comparison(self, table):
        """Get INF values grouped by MSA status"""
        query = f"""
        SELECT RNA_INF, af3_rna_MSA 
        FROM {table}
        WHERE RNA_INF IS NOT NULL AND af3_rna_MSA IS NOT NULL
        """
        rows = self.execute_query(query)
        
        inf_values = {
            '+MSA': {'all': [], 'wc': [], 'nonwc': [], 'stack': []},
            '-MSA': {'all': [], 'wc': [], 'nonwc': [], 'stack': []}
        }
        
        for inf, msa in rows:
            if inf and msa:
                try:
                    # Process INF values
                    if isinstance(inf, str):
                        inf_list = eval(inf)
                        if isinstance(inf_list, list) and len(inf_list) == 4:
                            # Only include valid INF values (-1 indicates invalid/missing)
                            inf_values_filtered = [v for v in inf_list if v != -1]
                            if len(inf_values_filtered) != len(inf_list):
                                continue  # Skip if any value was -1
                                
                            # Process MSA values
                            if isinstance(msa, str):
                                msa_values = eval(msa)
                                # Classify as -MSA if any value is 0, otherwise +MSA
                                is_msa = 0 if (isinstance(msa_values, list) and 0 in msa_values) else 1
                            else:
                                is_msa = int(msa)
                            
                            # Add to appropriate lists
                            msa_key = '+MSA' if is_msa == 1 else '-MSA'
                            inf_values[msa_key]['all'].append(float(inf_list[0]))
                            inf_values[msa_key]['wc'].append(float(inf_list[1]))
                            inf_values[msa_key]['nonwc'].append(float(inf_list[2]))
                            inf_values[msa_key]['stack'].append(float(inf_list[3]))
                            
                except (ValueError, SyntaxError) as e:
                    print(f"Error processing values - INF: {inf}, MSA: {msa}")
                    continue
        
        return inf_values

    def get_msa_di_comparison(self, table):
        """Get DI values grouped by MSA status"""
        query = f"""
        SELECT RNA_DI, af3_rna_MSA 
        FROM {table}
        WHERE RNA_DI IS NOT NULL AND af3_rna_MSA IS NOT NULL
        """
        rows = self.execute_query(query)
        
        di_values = {'+MSA': [], '-MSA': []}
        for di, msa in rows:
            if di and msa:
                try:
                    # Convert DI value to float
                    if isinstance(di, str):
                        di = float(di)
                    
                    # Process MSA values
                    if isinstance(msa, str):
                        msa_values = eval(msa)
                        # Classify as -MSA if any value is 0, otherwise +MSA
                        is_msa = 0 if (isinstance(msa_values, list) and 0 in msa_values) else 1
                    else:
                        is_msa = int(msa)
                    
                    # Add to appropriate list
                    if is_msa == 1:
                        di_values['+MSA'].append(float(di))
                    else:
                        di_values['-MSA'].append(float(di))
                except (ValueError, SyntaxError) as e:
                    print(f"Error processing values - DI: {di}, MSA: {msa}")
                    continue
        
        return di_values['+MSA'], di_values['-MSA']

class DNAMetrics(DatabaseMethods, PlotCreator):
    def __init__(self, id: str = "", all_rmsd: str = "", dna_rmsd: str = "", dna_lddt: float = None,
                 dna_tm_score: float = None):
        DatabaseMethods.__init__(self)
        PlotCreator.__init__(self, 'evaluation_metrics', msa_option, single_chain_only)
        self.id = id
        self.all_rmsd = all_rmsd
        self.dna_rmsd = dna_rmsd
        self.dna_lddt = dna_lddt
        self.dna_tm_score = dna_tm_score

    def load_from_db(self, select_table):
        columns = ["Complex_RMSD", "DNA_RMSD", "DNA_LDDT", "DNA_TM"]
        condition = "exp_db_id = ?"
        params = (self.id,)

        data = self.get_table_columns(select_table, columns, condition, params)

        if data:
            (
                self.all_rmsd,
                self.dna_rmsd,
                self.dna_lddt,
                self.dna_tm_score
            ) = data

    def get_dna_rmsd_values(self, select_table):
        rows = self.get_table_columns(select_table, ["DNA_RMSD"], condition=None)
        dna_rmsd = []
        for row in rows:
            if row[0]:
                try:
                    value_list = ast.literal_eval(row[0])
                    if isinstance(value_list, list) and len(value_list) > 0:
                        dna_rmsd.append(value_list[0])  # Take the first value
                except (ValueError, SyntaxError):
                    continue
            else:
                continue
        return dna_rmsd

    def get_dna_lddt_values(self, select_table):
        rows = self.get_table_columns(select_table, ["DNA_LDDT"], condition=None)
        lddt_values = [row[0] for row in rows if row[0] not in (None, '')]
        return lddt_values

    def get_dna_tm_values(self, select_table):
        rows = self.get_table_columns(select_table, ["DNA_TM"], condition=None)
        dna_tm_score = [row[0] for row in rows if row[0] not in (None, '')]
        return dna_tm_score

    def get_x_and_length(self, x, select_table):
        condition = f"{x} IS NOT NULL AND {x} != '' AND RNALength IS NOT NULL"
        data = self.get_table_columns(select_table, ["exp_db_id", x, "DNALength"], condition)
        return max_length(data)

    def get_rmsd_and_length(self, select_table):
        condition = "DNA_RMSD IS NOT NULL AND DNA_RMSD != '' AND DNALength IS NOT NULL"
        rows = self.get_table_columns_all(select_table, ["exp_db_id", "DNA_RMSD", "DNALength"], condition)
        processed_data = []
        for row in rows:
            exp_db_id = row[0]
            dna_rmsd_str = row[1]
            dna_length = row[2]
            try:
                dna_rmsd = float(dna_rmsd_str.strip('[]').split(',')[0]) #pymol: 0, tmscore:1
                processed_data.append((exp_db_id, dna_rmsd, dna_length))
            except (ValueError, IndexError):
                continue
        return max_length(processed_data)

def plot_metrics(table, msa_option=None, single_chain_only=False):
    """Plot evaluation metrics with filtering options"""
    # Create plot creator first
    plot_creator = PlotCreator('evaluation_metrics', msa_option, single_chain_only)

    # Initialize metrics classes with plot_creator
    rna_metrics = RNAMetrics(msa_option=msa_option, single_chain_only=single_chain_only)

    if table == "pred_protein_rna" or table == "pred_protein_rna_dna" or table == "pred_protein_dna":
        protein_metrics = ProteinMetrics(msa_option=msa_option, single_chain_only=single_chain_only)
        name = "Protein"
        prot_rmsd_values = protein_metrics.get_protein_rmsd_values(table)
        print("!!!protein_rmsd_values", prot_rmsd_values)
        lddt_values = protein_metrics.get_protein_lddt_values(table)
        print("protein lddt_values", lddt_values)
        tm_values = protein_metrics.get_protein_tm_values(table)
        print("protein tm_values", tm_values)
        data_values = lddt_values, tm_values
        plot_creator.get_violin_plot(table,
                                     data_values=[lddt_values, tm_values],
                                     labels=[f'{name} LDDT', f'{name} TM'],
                                     name = "Protein_LDDTtm"
                                    )

    if table == "pred_protein_rna" or table == "pred_rna_rna":
        name = "RNA"
        rna_rmsd_values = rna_metrics.get_rna_rmsd_values(table)
        print("rna_rmsd_values", rna_rmsd_values)
        lddt_values = rna_metrics.get_rna_lddt_values(table)
        print("rna lddt_values", lddt_values)
        tm_values = rna_metrics.get_rna_tm_values(table)
        print("rna tm_values", tm_values)
        plot_creator.get_violin_plot(table,
                                     data_values=[lddt_values, tm_values],
                                     labels=[f'{name} LDDT', f'{name} TM'],
                                     name="RNA_LDDTtm"
                                     )
        if table == "pred_rna_rna": #protein_metrics.get_protein_rmsd_values(table) == []:
            plot_creator.get_violin_plot(table,
                                         data_values=[rna_rmsd_values],
                                         labels=[f'{name} RMSD'],
                                         name="RNA_RMSD"
                                         )
        if table == "pred_protein_rna" and protein_metrics.get_protein_rmsd_values(table) != []:
            plot_creator.get_violin_plot(table,
                                         data_values=[prot_rmsd_values, rna_rmsd_values],
                                         labels=["Protein RMSD","RNA RMSD"],
                                         name="Protein_RNA_RMSD"
                                         )

        all_inf, wc_inf, nonWC_inf, stacking_inf = rna_metrics.get_rna_inf_values(table)
        all_inf = [value for value in all_inf if value != -1]
        wc_inf = [value for value in wc_inf if value != -1]
        nonWC_inf = [value for value in nonWC_inf if value != -1]
        stacking_inf = [value for value in stacking_inf if value != -1]
        plot_creator.get_boxplot(table,
                                 data_values=[all_inf, wc_inf, nonWC_inf, stacking_inf],
                                 labels=["INF", "INF-WC", "INF-nonWC", "INF-STACK"],
                                 name="RNA_INFs"
                                 )
        rna_di = rna_metrics.get_rna_di_values(table)
        print("all_inf, wc_inf, nonWC_inf, stacking_inf", all_inf, wc_inf, nonWC_inf, stacking_inf)
        print("rna_di", rna_di)
        di_length_data = rna_metrics.get_di_and_length(table)
        plot_creator.get_scatterplot(table,
                                     di_length_data,
                                     xAxis_label = "RNA Sequence Length",
                                     yAxis_label = "DI",
                                     name="DI_RNAlength")

        lddt_length_data = rna_metrics.get_x_and_length("RNA_LDDT", table)
        # get_scatterplot(lddt_length_data, "LDDT", table, "RNA Sequence Length")
        plot_creator.get_scatterplot(table,
                                     lddt_length_data,
                                     xAxis_label="RNA Sequence Length",
                                     yAxis_label="LDDT",
                                     name="LDDT_RNAlength")
        rmsd_length_data = rna_metrics.get_rmsd_and_length(table)
        plot_creator.get_scatterplot(table,
                                     rmsd_length_data,
                                     xAxis_label="RNA Sequence Length",
                                     yAxis_label="RMSD",
                                     name="RMSD_RNAlength")
        tm_length_data = rna_metrics.get_x_and_length("RNA_TM", table)
        plot_creator.get_scatterplot(table,
                                     tm_length_data,
                                     xAxis_label="RNA Sequence Length",
                                     yAxis_label="TM-score",
                                     name="TM_RNAlength")

        exp_ids, extra_vals, all_inf, wc_inf, nonWC_inf, stacking_inf = rna_metrics.get_rna_inf_values_with_extra(
            table, "RNA_LDDT")
        plot_creator.get_scatterplot(table,
                                     all_inf,
                                     xAxis_label="INF",
                                     yAxis_label="RNA LDDT",
                                     yAxis_score = extra_vals,
                                     name="RNA-LDDT_INF")
        plot_creator.get_scatterplot(table,
                                     wc_inf,
                                     xAxis_label="INF-WC",
                                     yAxis_label="RNA LDDT",
                                     yAxis_score=extra_vals,
                                     name="RNA-LDDT_INF-WC")
        plot_creator.get_scatterplot(table,
                                     nonWC_inf,
                                     xAxis_label="INF-nonWC",
                                     yAxis_label="RNA LDDT",
                                     yAxis_score=extra_vals,
                                     name="RNA-LDDT_INF-nonWC")
        plot_creator.get_scatterplot(table,
                                     stacking_inf,
                                     xAxis_label="INF-STACK",
                                     yAxis_label="RNA LDDT",
                                     yAxis_score=extra_vals,
                                     name="RNA-LDDT_INF-STACK")
        plot_creator.get_scatterplot(table,
                                     stacking_inf,
                                     xAxis_label="INF-STACK",
                                     yAxis_label="INF-nonWC",
                                     yAxis_score=nonWC_inf,
                                     name="RNA-INF-nonWC_stack")
        exp_ids, extra_vals, all_inf, wc_inf, nonWC_inf, stacking_inf = rna_metrics.get_rna_inf_values_with_extra(
            table, "RNA_TM")
        plot_creator.get_scatterplot(table,
                                     all_inf,
                                     xAxis_label="INF",
                                     yAxis_label="RNA TM",
                                     yAxis_score=extra_vals,
                                     name="RNA-TM_INF")
        plot_creator.get_scatterplot(table,
                                     nonWC_inf,
                                     xAxis_label="INF-nonWC",
                                     yAxis_label="RNA TM",
                                     yAxis_score=extra_vals,
                                     name="RNA-TM_INF-nonWC")
        plot_creator.get_scatterplot(table,
                                     stacking_inf,
                                     xAxis_label="INF-STACK",
                                     yAxis_label="RNA TM",
                                     yAxis_score=extra_vals,
                                     name="RNA-TM_INF-STACK")

        gc_lddt_content = rna_metrics.get_x_and_y("RNA_GC_Content", "RNA_LDDT", table)
        plot_creator.get_scatterplot(table,
                                     gc_lddt_content,
                                     xAxis_label="RNA GC-Content",
                                     yAxis_label="RNA LDDT",
                                     name="RNA-GC_LDDT")

        conditionGC = "RNA_GC_Content >= 0.5"
        gc_pTM_content = rna_metrics.get_x_and_y_condition("af3_rna_ipTM", "RNA_TM", conditionGC, table)
        plot_creator.get_scatterplot(table,
                                     gc_pTM_content,
                                     xAxis_label="RNA pTM",
                                     yAxis_label="RNA TM",
                                     name="RNA-GC>05_pTM_TM")

        au_lddt_content = rna_metrics.get_x_and_y("RNA_AU_Content", "RNA_LDDT", table)
        plot_creator.get_scatterplot(table,
                                     au_lddt_content,
                                     xAxis_label="RNA AU-Content",
                                     yAxis_label="RNA LDDT",
                                     name="RNA-AU_LDDT")

        gc_tm_content = rna_metrics.get_x_and_y("RNA_GC_Content", "RNA_TM", table)
        plot_creator.get_scatterplot(table,
                                     gc_tm_content,
                                     xAxis_label="RNA GC-Content",
                                     yAxis_label="RNA TM",
                                     name="RNA-GC_TM")

        au_tm_content = rna_metrics.get_x_and_y("RNA_AU_Content", "RNA_TM", table)
        plot_creator.get_scatterplot(table,
                                     au_tm_content,
                                     xAxis_label="RNA AU-Content",
                                     yAxis_label="RNA TM",
                                     name="RNA-AU_TM")

        # Add MSA comparison violin plots for all metrics
        metrics_to_compare = [
            ('RNA_RMSD', 'RMSD [Å]'),
            ('RNA_TM', 'TM-score'),
            ('RNA_LDDT', 'LDDT')
        ]

        for metric_col, ylabel in metrics_to_compare:
            msa_values, no_msa_values = rna_metrics.get_msa_metric_comparison(table, metric_col)
            if msa_values and no_msa_values:  # Only create plot if we have data for both cases
                plot_creator.get_violin_plot(
                    table_source='rna_metrics',
                    data_values=[msa_values, no_msa_values],
                    labels=['+MSA', '-MSA'],
                    name=f'{metric_col}_MSA_comparison',
                    yAxis_label=ylabel
                )

        # Add MSA comparison boxplots for INF metrics
        inf_values = rna_metrics.get_msa_inf_comparison(table)
        if inf_values['+MSA']['all'] and inf_values['-MSA']['all']:  # Check if we have data
            # Create boxplot for each INF type
            inf_types = ['all', 'wc', 'nonwc', 'stack']
            labels = ['INF', 'INF-WC', 'INF-nonWC', 'INF-STACK']

            for inf_type, label in zip(inf_types, labels):
                plot_creator.get_boxplot(
                    table_source='rna_metrics',
                    data_values=[inf_values['+MSA'][inf_type], inf_values['-MSA'][inf_type]],
                    labels=['+MSA', '-MSA'],
                    name=f'{label}_MSA_comparison',
                    yAxis_label=label
                )

        # Add MSA comparison for DI
        msa_di, no_msa_di = rna_metrics.get_msa_di_comparison(table)
        if msa_di and no_msa_di:  # Only create plot if we have data for both cases
            plot_creator.get_boxplot(
                table_source='rna_metrics',
                data_values=[msa_di, no_msa_di],
                labels=['+MSA', '-MSA'],
                name='DI_MSA_comparison',
                yAxis_label='Deformation Index'
            )

    if table == "pred_protein_dna" or table == "pred_protein_rna_dna":
        dna_metrics = DNAMetrics(msa_option=msa_option, single_chain_only=single_chain_only)
        name = "DNA"
        dna_rmsd_values = dna_metrics.get_dna_rmsd_values(table)
        print("dna_rmsd_values", rna_rmsd_values)
        lddt_values = dna_metrics.get_dna_lddt_values(table)
        print("dna lddt_values", lddt_values)
        tm_values = dna_metrics.get_dna_tm_values(table)
        print("dna tm_values", tm_values)
        plot_creator.get_violin_plot(table,
                                     data_values=[lddt_values, tm_values],
                                     labels=[f'{name} LDDT', f'{name} TM'],
                                     name="DNA_LDDTtm"
                                     )
        if rna_metrics.get_rna_rmsd_values(table) != []:
            plot_creator.get_violin_plot(table,
                                         data_values=[prot_rmsd_values, rna_rmsd_values, dna_rmsd_values],
                                         labels=["Protein RMSD", "RNA RMSD", "DNA RMSD"],
                                         name="Protein-RNA-DNA_RMSD"
                                         )
        lddt_length_data = dna_metrics.get_x_and_length("DNA_LDDT", table)
        plot_creator.get_scatterplot(table,
                                     lddt_length_data,
                                     xAxis_label="DNA Sequence Length",
                                     yAxis_label="LDDT",
                                     name="LDDT_DNAlength")
        tm_length_data = dna_metrics.get_x_and_length("DNA_TM", table)
        plot_creator.get_scatterplot(table,
                                     tm_length_data,
                                     xAxis_label="DNA Sequence Length",
                                     yAxis_label="TM",
                                     name="TM_DNAlength")
        rmsd_length_data = dna_metrics.get_rmsd_and_length(table)
        plot_creator.get_scatterplot(table,
                                     rmsd_length_data,
                                     xAxis_label="DNA Sequence Length",
                                     yAxis_label="RMSD",
                                     name="RMSD_DNAlength")

    if table == "pred_protein_rna" or table == "pred_protein_dna" or table == "pred_protein_rna_dna":
        protein_metrics.close_connection()
    if table == "pred_rna_rna":
        rna_metrics.close_connection()

if __name__ == "__main__":
    if len(sys.argv) not in [2, 3, 4]:
        print("Usage: python evaluationMetrics.py [pred table] [single_chain|all] [+MSA|-MSA]")
        sys.exit(1)

    table_name = sys.argv[1]
    single_chain_only = sys.argv[2] == 'single_chain' if len(sys.argv) > 2 else False
    msa_option = sys.argv[3] if len(sys.argv) > 3 else None

    if msa_option and msa_option.upper() not in ['+MSA', '-MSA']:
        print("MSA option must be either +MSA or -MSA")
        sys.exit(1)

    plot_metrics(table_name, msa_option, single_chain_only)
