from TGA_FTIR_tools import Sample, Worklist, Baseline
from pytest import fixture
from pathlib import Path
import matplotlib.pyplot as plt

cali_sample_names = {
    "Netzsch": [f'ExpDat_250605_dd_Calciumoxalat_0{i}_spline' for i in range(2, 8)],
    "Otto": [f"dd_220722_Ca-oxalat_0{i}" for i in range(1,2)]
}

sample_names = {
    "Netzsch": ["250919_dd_PAC+H2O_02", "250919_dd_PAC+Permeat_03", cali_sample_names["Netzsch"][0]],
    "Otto": ["dd_220721_AS5000_01", "dd_220722_FU1_01", cali_sample_names["Otto"][0]]
}


@fixture(params = [(profile, samples) for profile, samples in cali_sample_names.items()], ids=cali_sample_names.keys())
def get_cali_worklist(request):
    profile, names = request.param
    yield Worklist(names, profile=profile)

@fixture(params = [(profile, samples) for profile, samples in sample_names.items()], ids=sample_names.keys())
def get_sample(request):
    profile, names = request.param
    #TODO other init methods
    for name in names:  
        yield Sample(name, profile=profile)

@fixture()
def get_baseline(request):
    yield Baseline("dd_220721_auftriebsblindwert_600-03", profile="Otto")

@fixture(params = [(profile, samples) for profile, samples in sample_names.items()], ids=sample_names.keys())
def get_worklist(request):
    profile, names = request.param
    def init_options():
        yield Worklist(names, profile=profile)

        #TODO add from_pickle, from_samplelog
        for fun in [lambda x, y: Sample(x, profile=y), lambda x, y: [Sample(x, profile=y)]]:
            init = fun(*request.param)
            yield Worklist(init)
    return init_options

@fixture(params = [(profile) for profile, samples in sample_names.items()])
def get_options(request, get_baseline):
    profile, samples = request.param
    ax = plt.subplots()
    #baseline = (get_baseline, ) #Baseline()
    yield {
        'ax':(ax, None),
        'baseline':(None, ),
        'corrs':({}, {},),
        'directory':(None, Path(".")),
        'key':(None, ),
        'name': tuple(samples),
        'plot':(None, ),
        'presets':(None, ),
        'profile': tuple(sample_names.keys()),
        'ref_mass_name':(None, ),
        'reference_name':("test", "dittmann2021c"),
        'T_max':(600, 1200),
        'T_max_tol':(None, ),
        'values':(None, )
    }