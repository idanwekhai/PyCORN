import numpy as np
import pytest

from pycorn import PcUni6
from pycorn.utils import import_xml_as_df


def test_import_utility():
    file_path = r"..\samples\sample.zip"
    reference_data = np.array([[1.72348633e+01, 9.91999893e+01, -1.64260855e-03, 9.71418340e-03, 0.00000000e+00],
                               [1.69039273e+01, 9.91999893e+01, 7.54344136e-01, 8.22319886e-01, 0.00000000e+00],
                               [1.70131853e+01, 9.91999893e+01, 1.12997290e-02, 2.78141668e-02, 0.00000000e+00],
                               [1.70404669e+01, 9.91999893e+01, -4.40397997e-02, -7.03156144e-03, 0.00000000e+00],
                               [1.70849351e+01, 9.91999893e+01, -3.79689479e-02, -1.14569797e-02, 0.00000000e+00]])

    data = import_xml_as_df(file_path, index=np.arange(-14, 1, 3))
    # ((data - data.min()) / (data.max() - data.min())).plot()
    assert data.values == pytest.approx(reference_data)


def test_pcuni6_interface():
    file_path = r"..\samples\sample.zip"
    reference_data = [(0.0, -0.0016426085494458675),
                      (0.0, -0.0016426085494458675),
                      (0.0, -0.0016426085494458675),
                      (0.0, -0.0016426085494458675),
                      (0.0, -0.0016426085494458675),
                      (0.0, -0.002494621090590954),
                      (0.0, -0.0035596368834376335),
                      (0.0, -0.004624652676284313),
                      (0.0, -0.00568966893479228),
                      (0.0, -0.006754684261977673)]
    reference_keys = ['CalibrationSettingData', 'Chrom.1.Xml', 'Chrom.1_10_True', 'Chrom.1_13_True', 'Chrom.1_14_True',
                      'Chrom.1_15_True', 'Chrom.1_16_True', 'Chrom.1_1_True', 'Chrom.1_20_True', 'Chrom.1_21_True',
                      'Chrom.1_23_True', 'Chrom.1_2_True', 'Chrom.1_38_True', 'Chrom.1_3_True', 'Chrom.1_4_True',
                      'Chrom.1_5_True', 'Chrom.1_6_True', 'Chrom.1_7_True', 'Chrom.1_89_True', 'Chrom.1_8_True',
                      'Chrom.1_90_True', 'Chrom.1_9_True', 'ColumnIndividualData', 'ColumnTypeData',
                      'EvaluationLog.xml', 'EvaluationProcedureData', 'InstrumentConfigurationData', 'Manifest.xml',
                      'MethodData', 'MethodDocumentationData', 'NextBufferPrepData', 'NextFracData', 'ReportFormatData',
                      'Result.xml', 'StrategyData', 'SystemData', 'SystemSettingData', 'VersionInformationData',
                      'Chrom.1.Xml_dict']

    xml_data = PcUni6(file_path)
    xml_data.load()
    assert reference_keys == list(xml_data.keys())
    xml_data.load_all_xml()
    assert reference_data == xml_data["Chrom.1"]["UV 1_280"]["data"][:10]
