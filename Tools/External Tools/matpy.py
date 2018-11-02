# matpy.py
def datenum_to_datetime(datenums):
    import numpy as np
    import pandas as pd
    """
    Convert Matlab datenum into Python datetime.
    :param datenum: Date in datenum format
    :return:        Datetime object corresponding to datenum.
    """
    
    dates=pd.to_datetime(np.array(datenums)-719529, unit='D')
   
    return dates.to_pydatetime()