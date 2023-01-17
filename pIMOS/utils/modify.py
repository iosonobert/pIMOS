import numpy as np
import datetime
"""

    *imosModifyTime* 
    
"""

def imosModifyTime(rr, add_time, comment=None):
    """
    Offset time by +add_time and log in the ds comments.
    """
    
    rr.ds = rr.ds.assign_coords({'time': rr.ds.time.values+add_time})  
    
    string = 'Time variable offset by adding {} with user comment "{}"'.format(str(add_time), comment)
    
    print(string)

    rr.add_comment(rr.default_user, string) # Comment the main dataset
    rr.add_comment(rr.default_user, string, data_var='time') # Comment the time variable too