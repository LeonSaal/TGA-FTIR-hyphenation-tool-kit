from TGA_FTIR_tools import Sample, Worklist, Baseline
from pytest import fixture

cali_sample_names = {
    "Netzsch": [f'ExpDat_250605_dd_Calciumoxalat_0{i}_spline' for i in range(2, 8)],
    "Otto": [f"dd_220722_Ca-oxalat_0{i}" for i in range(1,2)]
}

sample_names = {
    "Netzsch": ["250919_dd_PAC+H2O_02", "250919_dd_PAC+Permeat_03", cali_sample_names["Netzsch"][0]],
    "Otto": ["dd_220721_AS5000_01", "dd_220722_FU1_01", cali_sample_names["Otto"][0]]
}

sample_names.keys() = ["Netzsch", "Otto"]

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