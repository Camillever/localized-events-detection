
class Station():

    def __init__(self, name, position, channels, network, sampling_rate, sensor):
        self.name = name
        self.position = position
        self.channels = channels
        self.network = network
        self.sampling_rate = sampling_rate
        self.sensor = sensor

    def __repr__(self):
        return f"{self.name}({self.network})"

FRE = Station(
    'FRE',
    sampling_rate=100,
    position=(-21.201837, 55.695393, 1775),
    channels=('HHE', 'HHN', 'HHZ'),
    network='PF',
    sensor='GURALP CMG-3ESPC'
)

HIM = Station(
    'HIM',
    sampling_rate=100,
    position=(-21.211121, 55.720224, 1958),
    channels=('HHE', 'HHN', 'HHZ'),
    network='PF',
    sensor='GURALP CMG-3ESPCDE 100S'
)

NSR = Station(
    'NSR',
    sampling_rate=100,
    position=(-21.206879, 55.723376, 2035),
    channels=('EHZ'),
    network='PF',
    sensor='MARCK PRODUCT L4C 1Hz'
)

PER = Station(
    'PER',
    sampling_rate=100,
    position=(-21.185722, 55.686033, 2085),
    channels=('EHZ'),
    network='PF',
    sensor='MARCK PRODUCT L4C 1Hz'
)

TTR = Station(
    'TTR',
    sampling_rate=100,
    position=(-21.188532, 55.779374, 856),
    channels=('EHZ'),
    network='PF',
    sensor='MARCK PRODUCT L4C 1Hz'
)

RER = Station(
    'RER',
    sampling_rate=20,
    position=(-21.1712, 55.73986, 834),
    channels=('BHN', 'BHN', 'BHZ'),
    network='G',
    sensor='GEOTECH SL-210'
)

STATIONS = [FRE, HIM, NSR, PER, TTR, RER]
STATIONS_NETWORKS = {}
for station in STATIONS:
    try:
        STATIONS_NETWORKS[station.network][station.name] = station
    except KeyError:
        STATIONS_NETWORKS[station.network] = {station.name: station}

