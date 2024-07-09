import simulate
gaussiancutoff = 0.5
# Currently must be "tophat" or "fred" or "gaussian" or "wilma" or "ered" or 'halfgaussian1' 'parabolic' or 'choppedgaussian'
simulate.read_ini_file()
simulate.initialise('choppedgaussian')
config = simulate.get_configuration()


