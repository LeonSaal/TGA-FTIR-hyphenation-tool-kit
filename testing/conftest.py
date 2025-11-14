from TGA_FTIR_tools import Sample, Worklist, Baseline
from pytest import fixture

cali_sample_names = [['ExpDat_250605_dd_Calciumoxalat_02_spline',
 'ExpDat_250605_dd_Calciumoxalat_03_spline',
 'ExpDat_250605_dd_Calciumoxalat_04_spline',
 'ExpDat_250605_dd_Calciumoxalat_05_spline',
 'ExpDat_250605_dd_Calciumoxalat_06_spline',
 'ExpDat_250605_dd_Calciumoxalat_07_spline'], ["dd_220722_Ca-oxalat_02", "dd_220722_Ca-oxalat_01"]]
sample_names = [cali_sample_names[0][0], "dd_220721_AS5000_01"]
profiles = ["Netzsch", "Otto"]

@fixture(params = [(samples, profile) for samples, profile in zip(cali_sample_names, profiles)], ids=profiles)
def get_cali_worklist(request):
    names, profile = request.param
    yield Worklist(names, profiles=profile)

@fixture(params = [(name, profile) for name, profile in zip(sample_names, profiles)], ids=profiles)
def get_sample(request):
    name, profile = request.param
    yield Sample(name, profile=profile)

@fixture()
def get_baseline(request):
    yield Baseline("dd_220721_auftriebsblindwert_600-03", profile="Otto")

@fixture(params = ((samplename, profile) for samplename, profile in zip(sample_names*3, profiles*3)), ids=profiles*3)
def get_worklist(request):
    def init_options():
        yield Worklist(request.param[0], profiles = request.param[1])
        for fun in [lambda x, y: Sample(x, profile=y), lambda x, y: [Sample(x, profile=y)]]:
            init = fun(*request.param)
            yield Worklist(init)
    return init_options