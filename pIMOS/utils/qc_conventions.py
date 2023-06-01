
class base():
    """
    Base class for QAQC conventions. 
    """

    convention = 'base'
    convention_version = 'a'
    description = 'Good data: 0; Bad data: positive values; Under investigation: negative values'
    data_format = 'numeric'

    def to_dict(self):
        """
        Dictionary of attributes for a QAQC variable.
        """

        dict_out = {}
        dict_out['convention'] = self.convention
        dict_out['convention_version'] = self.convention_version
        dict_out['data_format'] = self.data_format
        dict_out['description'] = self.description

        dict_out['is_qc_flag'] = 1
        dict_out['long_name'] = 'qc_flag'
        dict_out['cf_compliant'] = 0
        dict_out['associated_data_variables'] = '' # The program must grow this list. 

        return dict_out

    def run_compliance(self, da):
        """
        Writes the metadata for the QC conventions to the attributes of an xarray dataarray and checks that the dataarray is compliant with the conventions.

        """

        print('WARNING, COMPLIANCE CHECK NOT YET IMPLEMENTED.')

        da.attrs.update(self.to_dict())

