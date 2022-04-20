# Make cross section plots for all casts
from CastAway import *

# Build casts class
casts = read_CastAway_csv('./data/surveys/CTD/casts',
                          './data/surveys/CTD/meta/cast_meta.csv')

# Get survey/transect ids
sections = casts.return_section_names()
max_depth = round(casts.max_depth())
out_dir = "./output/xsections"


# 2018_S1_T2

fig_S1 = casts.plot_xsection('2018_S1', 'T2',
    variable="temperature", set_crossshore_dist=15, vmin=12, vmax=21,
    acoustics=False, set_depth=80)[0]

fig_S2 = casts.plot_xsection('2018_S2', 'T2',
    variable="temperature", set_crossshore_dist=15, vmin=12, vmax=21,
    acoustics=False, set_depth=80)[0]

fig_S1.set_size_inches(15, 4)
fig_S2.set_size_inches(15, 4)

fig_S1.savefig("./figures/xsection/2018_S1_T1.png", dpi=300)
fig_S2.savefig("./figures/xsection/2018_S2_T1.png", dpi=300)