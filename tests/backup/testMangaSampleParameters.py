from hooloovookit import DatabaseConnection
a = DatabaseConnection.DatabaseConnection(user='albireo', name='apo_platedb')

from mangadb import SampleModelClasses as sModel

session = a.Session()
session.begin(subtransactions=True)

# nn = session.query(sModel.Parameter).get(1)
# nn.parameter = 'gender'
pT = session.query(sModel.ParentTargetToValue).all()
print(pT)
pT = session.query(sModel.ParentTarget).all()
print(pT[0].data)

# pT[0].data['height'] = '1.71'
# del pT[0].data['height']

pT = session.query(sModel.ParentTargetToValue).all()
print(pT)

# pT = session.query(sModel.Parameter).all()
# for pp in pT:
#     print(pp)

# newValue = sModel.Value(value='blue', parameter='eyes')
# session.add(newValue)
session.commit()

print(session.query(sModel.Value).all())
