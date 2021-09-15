# Make cross section plots for all casts
from CastAway import *

# Build casts class
casts = read_CastAway_csv('./data/surveys/CTD/casts',
                          './data/surveys/CTD/meta/cast_meta.csv')

# Get survey/transect ids
sections = casts.return_section_names()
max_depth = round(casts.max_depth())
out_dir = "./output/xsections"

# # make plots
# for var in ['temperature', 'salinity']:
#     vrange = casts.variable_range(var)
#     for i, row in sections.iterrows():
#         # make plots
#         figs = casts.plot_xsection(row.survey, row.transect,
#                                     variable=var,
#                                     set_crossshore_dist=15,
#                                     vmin=vrange[0], vmax=vrange[1],
#                                     set_depth=max_depth,
#                                     acoustics=True)
#         # save plots
#         for j in range(0, len(figs)):
#             figs[j].savefig(out_dir+"/"+row.survey+"_"+row.transect+"_"+var+"_run"+str(j), dpi=300)

# Return swarm enviromental data
agg = casts.swarm_environment_data()
agg = agg.reset_index(drop=True)
agg.to_csv('./data/surveys/acoustics/agg_CTD_interporlations.csv', index=False)