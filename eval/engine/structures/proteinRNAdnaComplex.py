from database.startConfig import StartConfig
from database.databaseMethods import DatabaseMethods

class RNAProteinDNAcomplex:
    config = StartConfig()
    pred_table_name = config.pred_table
    exp_table_name = config.ref_table
    def __init__(self, id: str, is_pred: bool, file_name: str = None,
                 contact_list: str = "",
                 interface_atom_ids: str = "", interface_rna_atom_ids: str = "", dnaProt_interface_atom_ids: str = "",
                 interface_dna_atom_ids: str = ""):
        self.id = id
        self.is_pred = is_pred
        self.file_name = file_name if is_pred else None
        self.interface_metrics = []
        self.contact_list = contact_list
        self.interface_atom_ids = interface_atom_ids
        self.interface_rna_atom_ids = interface_rna_atom_ids
        self.dnaProt_interface_atom_ids = dnaProt_interface_atom_ids
        self.interface_dna_atom_ids = interface_dna_atom_ids

    @classmethod
    def choose_interface_atom_ids(cls, row):
        if cls.pred_table_name in ['pred_protein_rna', 'pred_protein_rna_dna']:
            return row['rna_prot_interface_atom_ids']
        elif cls.pred_table_name == 'pred_protein_dna':
            return row['dna_prot_interface_atom_ids']
        elif cls.pred_table_name == 'pred_rna_rna':
            return row['rna_rna_interface_atom_ids']
        else:
            raise ValueError("ERROR: Invalid table name for choosing interface_atom_ids")

    @classmethod
    def choose_interface_na_atom_ids(cls, row):
        if cls.pred_table_name in ['pred_protein_rna', 'pred_protein_rna_dna', 'pred_rna_rna']:
            return row['interface_rna_atom_ids']
        elif cls.pred_table_name == 'pred_protein_dna':
            return row['interface_dna_atom_ids']
        else:
            raise ValueError("ERROR: Invalid table name for choosing interface_na_atom_ids")

    @classmethod
    def get_proteinRNAdnaComplex_from_db_predTable(cls, id: str, file_name=None):
        database_methods = DatabaseMethods()
        id = id.upper()
        # if file_name == None:
        #     table_name = cls.exp_table_name
        # else:
        table_name = cls.pred_table_name
        row = database_methods.get_table_row(table_name = table_name,
                                             condition = f"exp_db_id = '{id}' AND FileName = '{file_name}'")
        if row is None:
            raise ValueError(f"No row was selected in table {table_name} for id {id} and file name {file_name}.")

        # if table_name == 'exp_protein_rna':
        #     proteinRNAdnaComplex = RNAProteinDNAcomplex(
        #         id=row['PDBId'],
        #         is_pred=True,
        #         file_name=row['FileName'],
        #         contact_list=row['ContactList'],
        #         interface_atom_ids=cls.choose_interface_atom_ids(row),
        #         interface_rna_atom_ids=cls.choose_interface_na_atom_ids(row)
        #     )
        if table_name != 'pred_protein_rna_dna':
            proteinRNAdnaComplex = RNAProteinDNAcomplex(
                id=row['exp_db_id'],
                is_pred=True,
                file_name=row['FileName'],
                contact_list=row['ContactList'],
                interface_atom_ids=cls.choose_interface_atom_ids(row),
                interface_rna_atom_ids=cls.choose_interface_na_atom_ids(row)
            )
        else:
            proteinRNAdnaComplex = RNAProteinDNAcomplex(
                id=row['exp_db_id'],
                is_pred=True,
                file_name=row['FileName'],
                contact_list=row['ContactList'],
                interface_atom_ids=cls.choose_interface_atom_ids(row),
                interface_rna_atom_ids=cls.choose_interface_na_atom_ids(row),
                dnaProt_interface_atom_ids=row['dna_prot_interface_atom_ids'],
                interface_dna_atom_ids=row['interface_dna_atom_ids']
            )
        # database_methods.close_connection()
        return proteinRNAdnaComplex

    def get_af3Metrics(self):
        return self.ptm, self.iptm, self.fraction_disordered

    def get_interfaceAtom_ids(self):
        return self.interface_atom_ids

    def get_RNA_interfaceAtom_ids(self):
        return self.interface_rna_atom_ids

    def get_dnaProt_interfaceAtom_ids(self):
        return self.dnaProt_interface_atom_ids

    def get_DNA_interfaceAtom_ids(self):
        return self.interface_dna_atom_ids

    def get_contact_list(self):
        return self.contact_list


# rnaProtein_instance = RNAProteinComplex.get_proteinRNAcomplex_from_db(id="2x1f", is_pred=True, file_name='fold_2x1f_s1040878328_model_0.cif')
# print("rnaProtein_instance", rnaProtein_instance.get_interfaceAtom_ids())
