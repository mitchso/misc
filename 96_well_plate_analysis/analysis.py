# TODO: Improve coloration of graph
# TODO:

import re
import string
import pandas
from sys import argv

import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib import rcParams
from scipy.stats import variation, tstd


class Sample(object):
    def __init__(self):
        self.name = "not assigned"
        self.wells = sorted([])
        self.values = []
        self.treatment = "not assigned"
        self.concentration = 0.0

    def average_od(self):
        return sum(self.values) / len(self.values)

    def cv(self):
        if len(self.values) > 1:
            return variation(self.values)*100
            # return tstd(self.values) * 100
        else:
            return 0

    def viability(self, no_tmt_od: float):
        return self.average_od() * 100 / no_tmt_od


class Condition(object):
    def __init__(self):
        self.name = "not assigned"
        self.samples = []

    def is_control(self):
        if re.search('cells', self.name):
            return True
        else:
            return False

    def max_concentration(self):
        return max([x.concentration for x in self.samples])

    def min_concentration(self):
        return min([x.concentration for x in self.samples])

    def sort_samples(self):
        return sorted(self.samples, key=lambda x: x.concentration)

    def color(self):
        return None  # TODO: Fix this if needed

    def style_parameters(self):
        return {
            'linestyle': 'dashed' if 'pep' in self.name.lower() else 'solid',
            'label': self.name,
            'capsize': 5,
            'marker': '.',
            'color': self.color(),
            'linewidth': 1.75
        }


class Experiment(object):
    def __init__(self):
        self.name = "not assigned"
        self.conditions = []
        self.controls = []

    def max_concentration(self):
        return max([x.max_concentration() for x in self.conditions])

    def min_concentration(self):
        return min([x.max_concentration() for x in self.conditions])

    def sort_conditions(self):
        return sorted(self.conditions, key=lambda x: x.name)


xtt_file = argv[1]
info_file = argv[2]
out_file = argv[3]
experiment = Experiment()

# Read info from info_file into memory
info_table = pandas.read_csv(info_file, sep='\t')
treatment_table = info_table.loc[info_table.ne(0).all(axis=1)]
control_table = info_table.loc[info_table['Treatment'].isin(['cells only (lysis)', 'cells only'])]
final_table = pandas.concat([treatment_table, control_table])
final_table['Well'] = final_table['Row']+final_table['Column'].map(str)

treatments = list(set([row[1]['Treatment'] for row in final_table.iterrows()]))
concentrations = list(set([row[1]['Concentration'] for row in final_table.iterrows()]))


# Collect sample information
sample_dict = {}
for treatment in treatments:
    for concentration in concentrations:
        sample = Sample()
        for index, row in final_table.iterrows():
            if row['Treatment'] == treatment:
                if row['Concentration'] == concentration:
                    sample.name = treatment+", "+str(concentration)+" nM"
                    sample.treatment = treatment
                    sample.concentration = concentration
                    sample.wells.append(row['Well'])
                    sample_dict[sample.name] = sample


# Fix weird encoding from old plate reader
fixed_file = []
with open(xtt_file, 'rb') as xtt_file:
    for line in xtt_file.readlines():
        fixed_1 = ''.join(filter(lambda x: x in string.printable, str(line)))   # remove non-ascii characters
        fixed_2 = fixed_1.replace('b\'', '')    # remove from start of each line
        fixed_3 = fixed_2.replace('\\t', '\t')  # remove non-ascii tabs
        fixed_4 = fixed_3.replace('\\r\\n\'', '\n')  # remove from end of each line
        fixed_file.append(fixed_4)  # this should be actually fixed


# Collect info from input file
for key in sample_dict:
    for well in sample_dict[key].wells:
        for line in fixed_file:
            if re.search(well, line):
                if line.split(sep='\t')[3] == 'Outlier':
                    od_value = line.split(sep='\t')[2]
                else:
                    od_value = line.split(sep='\t')[3]
                sample_dict[key].values.append(float(od_value))

# Need this for making calculations
no_treatment_od = sample_dict['cells only, 0.0 nM'].average_od()

# Group values into conditions
for treatment in treatments:
    condition = Condition()
    condition.name = treatment
    condition.samples = [sample for sample in sample_dict.values() if sample.treatment is treatment]
    experiment.conditions.append(condition)

# Write output table
open(out_file, 'w').close()
with open(out_file, 'w') as out_file:
    out_file.write("\t".join(["Treatment", "Concentration (nM)", "Wells", "Raw Values", "Mean OD", "CV %", "Viability", "\n"]))
    for condition in experiment.sort_conditions():
        for sample in condition.sort_samples():
            out_file.write("\t".join(str(x) for x in ([sample.treatment,
                                                       sample.concentration,
                                                       ", ".join(str(x) for x in sample.wells),
                                                       ", ".join(str(x) for x in sample.values),
                                                       round(sample.average_od(), 3),
                                                       round(sample.cv(), 3),
                                                       round(sample.viability(no_treatment_od), 3),
                                                       "\n"])))

print()

# Add data to graph
y_offset = 1.0
for condition in experiment.sort_conditions():
    if not condition.is_control():
        x = [sample.concentration*y_offset for sample in condition.sort_samples()]
        y = [sample.viability(no_treatment_od) for sample in condition.sort_samples()]
        errors = [sample.cv() for sample in condition.sort_samples()]
        #plt.plot(x, y, **condition.style_parameters())
        plt.errorbar(x, y, yerr=errors, **condition.style_parameters())

        #y_offset += 0.015


# Stylize the graph
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

print(rcParams['font.sans-serif'])

plt.xscale('log')
fig, ax = plt.gcf(), plt.gca()  # Get current figure and axes instance

ax.xaxis.set_major_formatter(ScalarFormatter())  # Need this to make the axis numbers display in plain
ax.yaxis.set_major_formatter(ScalarFormatter())  # rather than scientific notation
[ax.spines[spine].set_visible(False) for spine in ax.spines  # Remove borders around the top and right side
 if ax.spines[spine].spine_type is 'right'
 or ax.spines[spine].spine_type is 'top']

plt.hlines(100, xmin=experiment.min_concentration(), xmax=experiment.max_concentration(), linestyles='dashed')
plt.hlines(0, xmin=experiment.min_concentration(), xmax=experiment.max_concentration(), linestyles='dashed')
# ax.spines[spine].set_visible(False)
plt.legend(frameon=False)  # Add legend

# Finally, display the graph
plt.show()
