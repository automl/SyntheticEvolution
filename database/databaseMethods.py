import sqlite3
from database.startConfig import StartConfig
import pandas as pd
# import traceback

class DatabaseMethods():
    def __init__(self, *args, **kwargs):
        config = StartConfig()
        self.db_path = config.db_path
        # self.db_conn = sqlite3.connect(self.db_path)
        # self.cursor = self.db_conn.cursor()
        self.connection = sqlite3.connect(self.db_path)
        self.connection.row_factory = sqlite3.Row  # Allow dictionary-style access
        self.cursor = self.connection.cursor()

    def close_connection(self):
        # if self.db_conn:
        #     self.db_conn.close()
        # traceback.print_stack()
        if self.connection:
            print("Closing database connection.")
            self.connection.close()

    # Gets data based on query and optional parameters
    def get_data(self, query, params=()):
        self.cursor.execute(query, params)
        return self.cursor.fetchone()

    # Updates rows
    def update_data(self, table_name, columns, values, condition):
        set_clause = ', '.join([f"{col} = ?" for col in columns])
        query = f"UPDATE {table_name} SET {set_clause} WHERE {condition}"
        self.cursor.execute(query, values)
        self.connection.commit()
        return self.cursor.rowcount

    # Inserts data into table_name
    def insert_data(self, table_name, columns, values):
        placeholders = ', '.join(['?'] * len(columns))
        query = f"INSERT INTO {table_name} ({', '.join(columns)}) VALUES ({placeholders})"
        self.cursor.execute(query, values)
        self.connection.commit()

    # Gets all specified rows with all columns from table_name matching the condition.
    def get_table_row(self, table_name, condition):
        query = f"SELECT * FROM {table_name} WHERE {condition}"
        return self.get_data(query)

    def get_table_columns_all(self, table_name: str, columns: list, condition: str, params: tuple = ()):
        query = f"SELECT {', '.join(columns)} FROM {table_name} WHERE {condition}"
        with sqlite3.connect(self.db_path) as connection:
            self.cursor.execute(query, params)
            return self.cursor.fetchall()
    #
    # def get_table_columns(self, table_name: str, columns: list, condition: str = None, params: tuple = ()):
    #     query = f"SELECT {', '.join(columns)} FROM {table_name}"
    #     if condition:
    #         query += f" WHERE {condition}"
    #     with sqlite3.connect(self.db_path) as connection:
    #         self.cursor.execute(query, params)
    #         return self.cursor.fetchall()

    def get_table_columns(self, table_name: str, columns: list, condition: str = None, params: tuple = ()):
        query = f"SELECT {', '.join(columns)} FROM {table_name}"
        if condition:
            query += f" WHERE {condition}"
            with sqlite3.connect(self.db_path) as connection:
                self.cursor.execute(query, params)
                return self.cursor.fetchone()
        with sqlite3.connect(self.db_path) as connection:
            self.cursor.execute(query, params)
            return self.cursor.fetchall()

    def execute_query(self, query, params=None):
        with sqlite3.connect(self.db_path) as conn:
            cursor = conn.cursor()
            if params:
                cursor.execute(query, params)
            else:
                cursor.execute(query)
            rows = cursor.fetchall()
        return rows

    # Updates or inserts new data
    def update_or_insert(self, table_name, columns, values, condition):
        if self.update_data(table_name, columns, values, condition) == 0:
            self.insert_data(table_name, columns, values)

    # Adds new column between 2 specified columns
    def add_column_between(self, table_name, new_column_name, new_column_type, column_a, column_b):
        # Get the table's current schema and columns
        self.cursor.execute(f"PRAGMA table_info({table_name})")
        columns = [col[1] for col in self.cursor.fetchall()]

        if column_a not in columns or column_b not in columns:
            raise ValueError(f"Columns {column_a} and {column_b} must both exist in the table.")

        # Determine the new column order with the inserted column
        new_columns = []
        for col in columns:
            new_columns.append(col)
            if col == column_a:
                new_columns.append(new_column_name)

        # Define the new column definitions
        column_defs = [f"{col} TEXT" if col != new_column_name else f"{new_column_name} {new_column_type}" for col in
                       new_columns]

        # Start a transaction to create a new table with the modified schema
        self.cursor.execute("BEGIN")
        try:
            temp_table_name = f"{table_name}_temp"
            self.cursor.execute(f"ALTER TABLE {table_name} RENAME TO {temp_table_name}")
            self.cursor.execute(f"CREATE TABLE {table_name} ({', '.join(column_defs)})")

            # Insert data from the old table into the new table, filling the new column with NULL
            insert_cols = ", ".join([f"{col}" for col in new_columns])
            select_cols = ", ".join([f"{col}" if col != new_column_name else "NULL" for col in new_columns])
            self.cursor.execute(f"INSERT INTO {table_name} ({insert_cols}) SELECT {select_cols} FROM {temp_table_name}")

            self.cursor.execute(f"DROP TABLE {temp_table_name}")
            self.connection.commit()

        except sqlite3.Error as e:
            self.connection.rollback()
            raise e

    def drop_all_tables(self):
        try:
            with sqlite3.connect(self.db_path) as conn:
                cursor = conn.cursor()
                cursor.execute("PRAGMA foreign_keys = OFF;")  # Temporarily disable foreign key checks
                cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
                tables = cursor.fetchall()

                for (table_name,) in tables:
                    cursor.execute(f"DROP TABLE IF EXISTS {table_name}")

                conn.commit()
        except sqlite3.Error as e:
            print("An error occurred while dropping tables:", e)

    # Finds common IDs using INTERSECT between two tables
    def get_intersect_values(self, table1: str, column1: str, table2: str, column2: str,
                          condition1: str = None, condition2: str = None):

        condition1 = f"WHERE {condition1}" if condition1 else ""
        condition2 = f"WHERE {condition2}" if condition2 else ""

        query = f"""
                SELECT {column1} FROM {table1} {condition1}
                INTERSECT
                SELECT {column2} FROM {table2} {condition2}
            """

        self.cursor.execute(query)
        common_values = [row[0] for row in self.cursor.fetchall()]
        return common_values

    def get_row(self, table_name, columns=None, condition=None):
        """
        Fetch a single row from the database
        Args:
            table_name: name of the table
            columns: list of column names or None for all columns
            condition: WHERE clause condition or None
        Returns: 
            Dictionary with column names as keys
        """
        try:
            cols = '*' if columns is None else ', '.join(columns)
            query = f"SELECT {cols} FROM {table_name}"
            if condition:
                query += f" WHERE {condition}"
            query += " LIMIT 1"

            self.cursor.execute(query)
            row = self.cursor.fetchone()
            
            if row:
                # If using *, get column names from cursor description
                if columns is None:
                    columns = [desc[0] for desc in self.cursor.description]
                return {columns[i]: row[i] for i in range(len(columns))}
            return None
            
        except sqlite3.Error as err:
            print(f"Error fetching data: {err}")
            return None

    def select(self, table, columns, condition=None):
        """
        Select multiple rows from the database
        Returns: list of tuples containing the requested columns
        """
        try:
            query = f"SELECT {', '.join(columns)} FROM {table}"
            if condition:
                query += f" WHERE {condition}"

            self.cursor.execute(query)
            return self.cursor.fetchall()
            
        except sqlite3.Error as err:
            print(f"Error selecting data: {err}")
            return None

    def filter_by_msa(self, df, msa_option):
        """Filter dataframe based on MSA option"""
        if msa_option not in ['+MSA', '-MSA']:
            return df

        def get_msa_value(msa_str):
            try:
                if pd.isna(msa_str):
                    return None
                if isinstance(msa_str, str):
                    values = eval(msa_str)
                    return 1 if any(v == 1 for v in values) else 0
                return msa_str
            except:
                return None

        df['msa_value'] = df['af3_rna_MSA'].apply(get_msa_value)
        msa_filter = df['msa_value'] == (1 if msa_option == '+MSA' else 0)
        return df[msa_filter]

    def filter_by_singleChain(self, exp_df, pred_df):
        """Filter dataframes to keep only single chain pairs"""
        print("\nFiltering for single chain pairs...")
        # Filter experimental data
        filtered_exp = exp_df[
            (exp_df['ProteinChainIDs'].str.len() == 1) &
            (exp_df['RNAChainIDs'].str.len() == 1)
        ]
        # Filter predicted data
        filtered_pred = pred_df[
            pred_df['exp_db_id'].isin(filtered_exp['PDBId'])
        ]
        print(f"Remaining entries after single chain filtering: {len(filtered_exp)}")
        return filtered_exp, filtered_pred