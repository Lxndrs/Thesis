import json
import glob

savefiles = glob.glob("*.save")

for savefile in savefiles:
    data = json.load(open(savefile))
    combined_data = {k: [] for k in data.keys()}
    varfit_pattern = savefile.replace('positionswill',
                                      'positionssamp*')
    combined_file = savefile.replace('positionswill', 'variations')
    varfitfiles = glob.glob("Multi-Fit/samp*/" + varfit_pattern)
    for varfitfile in varfitfiles:
        vardata = json.load(open(varfitfile))
        for k, v in vardata.items():
            combined_data[k].append(v)
    with open(combined_file, 'w') as f:
        json.dump(combined_data, f)

        
