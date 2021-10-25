import os

absolute = os.path.abspath(os.getcwd())

data = {
    "data": absolute + '/../data/',
    "disprot": absolute + '/../data/disprot/',
    "mobidb": absolute + '/../data/mobidb/',
    "merged": absolute + '/../data/merged/',
    "functions": absolute + '/../data/disprot/functions/',
    "analyze": absolute + '/../data/analyze/',
}
