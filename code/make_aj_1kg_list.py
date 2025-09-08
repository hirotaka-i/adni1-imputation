import pandas as pd
onekg = pd.read_csv('data/all_hg38_filtered_chrpos_pop.txt', sep=r'\s+')
# Population dictionary for grouping
pop_dict = {
    'AMR': ['MXL','CLM','PEL','PUR'],
    'EAS': ['JPT','CDX','CHB','CHS','KHV','CHD'],
    'EUR': ['TSI','IBS','GBR','CEU', 'AJ', 'FIN'],
    'SAS': ['PJL','ITU','STU','GIH','BEB'],
    'AFR': ['GWD','MSL','ESN','GWJ','YRI','LWK','GWF','GWW'],
    'AAC': ['ASW','ACB'],
}
onekg['#FID'] = 0
onekg['label'] = onekg['Population']
onekg['label_group'] = onekg['Population'].map({pop: grp for grp, pops in pop_dict.items() for pop in pops})
aj_eur = pd.read_csv('temp/AJ/pop_split_out/with_1kg_mah_EUR.list', sep=r'\s+', header=None, names=['#FID', 'IID'])
aj_eur['label'] = 'AJ'
aj_eur['label_group'] = 'EUR'
merged = pd.concat([onekg[['#FID', 'IID', 'label', 'label_group']], aj_eur], ignore_index=True)
merged_eur = merged[merged['label_group'] == 'EUR']
print(merged_eur['label'].value_counts())
merged_eur[['#FID', 'IID']].to_csv('temp/AJ/aj_1kg_eur.list', sep='\t', index=False, header=False)
merged_eur[['#FID', 'IID', 'label']].to_csv('temp/AJ/aj_1kg_eur_label.list', sep='\t', index=False, header=True)
print(f'Done: {len(merged_eur)} samples (AJ + 1kg EUR)')
print(f'saved to: temp/AJ/aj_1kg_eur.list')
print(f'EUR label saved to: temp/AJ/aj_1kg_eur_label.list')