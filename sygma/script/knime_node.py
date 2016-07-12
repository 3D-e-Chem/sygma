import sygma
import pandas as pd

phase1_cycles = options['phase1_cycles']
phase2_cycles = options['phase2_cycles']

scenario = sygma.Scenario([
    [sygma.ruleset['phase1'], phase1_cycles],
    [sygma.ruleset['phase2'], phase2_cycles]
])

parents = input_table[parents_column_name]

metabolites = []
for parent in parents:
    metabolites_network = scenario.run(parent)
    metabolites_network.calc_scores()
    metabolites_network.add_coordinates()
    metabolites += metabolites_network.to_list()

metabolites_df = pd.DataFrame(metabolites)
output_table = pd.merge(input_table,
                        metabolites_df,
                        left_on=parents_column_name,
                        right_on='parent')
