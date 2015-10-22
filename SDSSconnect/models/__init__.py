import warnings

# warnings.filterwarnings(
#     'ignore',
#     'This declarative base already contains a class with the same '
#     'class name and module name as sqlalchemy.ext.automap')
warnings.filterwarnings('ignore',
                        'Skipped unsupported reflection of expression-based')

from relationships import createRelationships
