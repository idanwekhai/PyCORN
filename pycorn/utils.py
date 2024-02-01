from pathlib import Path

import numpy as np
import pandas as pd

from pycorn import PcUni6


def get_series_from_data_dict(data_dictionary, target_key, data_key_list):
    try:
        # select the first injection as the injection timestamp
        inject_timestamp = data_dictionary[target_key]["Injection"]["data"][-1][0]
    except KeyError:
        inject_timestamp = 0

    data_series_list = []
    for data_key in data_key_list:
        data_array = np.array(data_dictionary[target_key][data_key]["data"]).astype(float)
        data_series = pd.Series(data=data_array[:, 1], index=data_array[:, 0])
        # remove duplicates
        data_series = data_series[~data_series.index.duplicated()]
        # offset by the infection_timestamp
        data_series.index -= inject_timestamp

        data_series_list.append(data_series)

    df = pd.concat(data_series_list, axis=1)
    df.columns = data_key_list
    return df


def import_xml_as_df(file_path: (str | Path), data_key_list: list = None, index: np.ndarray = None) -> pd.DataFrame:
    """
    Import the contents of a Unicorn Res/zip file into a pd.Dataframe

    Parameters
    ----------
    file_path : str or Path, Path to the res or zip file.
    data_key_list: list, optional, Keys to include in the DataFrame. Default: ["Cond", "UV", "Conc B"]
    index: np.ndarray, optional, Array of shape (1, ), to be used as index in the returned pd.DataFrame.
        Units are the same as the original data.

    Returns
    -------
    dataframe : pd.DataFrame

    """
    if data_key_list is None:
        data_key_list = ["Cond", "UV", "Conc B"]

    data_dictionary = PcUni6(file_path)
    data_dictionary.load_all_xml()

    target_key_list = [
        key for key in data_dictionary
        if "events" not in key
           and "Cond" in key
           and any(s.lower() in key.lower() for s in ["Tracer", "Injection", "Chrom", "Breakthrough", "Elution"])
           and any(["UV" in sub_key for sub_key in data_dictionary[key]])
    ]

    if any("breakthrough" in key.lower() for key in target_key_list):
        target_key_list = [key for key in target_key_list if "breakthrough" in key.lower()]

    if len(target_key_list) == 0:
        return None

    if "UV" not in data_dictionary[target_key_list[0]]:
        data_key_list.remove("UV")
        data_key_list.extend([sub_key for sub_key in data_dictionary[target_key_list[0]] if
                              ("UV" in sub_key and not "cell path" in sub_key)])

    series_list = [get_series_from_data_dict(data_dictionary, chrom, data_key_list) for chrom in target_key_list]

    if index is None:
        index = np.linspace(series_list[0].index.min(), series_list[0].index.max(), 100).round(3)

    # align and unify the index
    index = index[index < series_list[0].index.max()]

    if len(data_key_list) > 1:
        column_names = [(target_key, data_key) for target_key in target_key_list for data_key in data_key_list]
        column_names += [("temporary_insert", 0)]
    else:
        column_names = target_key_list
        column_names += ["temporary_insert"]

    # combine datapoints into one DataFrame and interpolate onto the unified index
    series_list.append(pd.Series(data=np.NaN, index=index))
    dataframe = pd.concat(series_list, axis=1, join="outer")
    dataframe.columns = column_names
    dataframe = dataframe.interpolate("index", )
    dataframe = dataframe.loc[index, column_names[:-1]]

    return dataframe
