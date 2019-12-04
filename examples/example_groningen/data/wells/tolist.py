# coding: utf-8
import pandas as pd
t = pd.read_csv(r'Onttrekkingen_MIPWA_Stationair.csv')
sqs = t.loc[:, [c for c in t.columns if c.startswith('SQ')]]
sqs.columns = [int(c.replace('SQ', '')) for c in sqs.columns]
sqs = sqs.stack()
t.index.names = 'id',
sqs.index.names = 'id', 'ilay'
sqs.name = 'q'
t2 = t.loc[:, [c for c in t.columns if not c.startswith('SQ')]].join(sqs, how='left')
t2.reset_index().to_csv('sq_list.csv')

