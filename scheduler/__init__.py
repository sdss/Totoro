
import warnings
from Totoro.exceptions import TotoroUserWarning

try:
    from observingPlan import ObservingPlan
    observingPlan = ObservingPlan()
except:
    warnings.warn('No observing plan found. Some functionality will be '
                  'limited.', TotoroUserWarning)
    observingPlan = None


from scheduler import BaseScheduler, Planner, Nightly, Plugger
