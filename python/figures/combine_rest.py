import pandas as pd

samples = pd.read_csv(snakemake.input['samples'])

no_movement = []
for index, sample in samples.iterrows():
    animal_id = sample['AnimalID']
    cell_name = sample['CellName']
    prefix = f'{snakemake.params["prefix"]}/{animal_id}/{cell_name}'
    dt = pd.read_csv(f'{prefix}/no_movement_final.csv')
    dt['AnimalID'] = animal_id
    dt['CellName'] = cell_name
    dt['SID'] = sample['SID']
    dt['Region'] = sample['Region']
    dt['EventID'] = [f'Event{x}' for x in range(len(dt))]
    no_movement.append(dt)

no_movement = pd.concat(no_movement, ignore_index = True)

no_movement.to_csv(snakemake.output['data'])
