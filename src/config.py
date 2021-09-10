import os

absolute = os.path.abspath(os.getcwd())

data = {
    "data": absolute + '/../data/',
    "disprot": absolute + '/../data/disprot/',
    "mobidb": absolute + '/../data/mobidb/',
    "merged": absolute + '/../data/merged/',
    "analyze": absolute + '/../data/analyze/',
}
