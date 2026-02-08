# Simulation scripts

These scripts use the [cells](https://github.com/yketa/cells) package.

## Epiboly simulations

We perform epiboly (stretching) simulations with `stretch.py` which takes an input open boundary configuration. We provide a standard simulation output file `disc.p` which generates such a configuration, and makes it circular and disordered using boundary tension and active Brownian driving.

```python
python stretch.py disc.p            # run stretching simulation
```

## Cutting/closing simulations

Starting from an initially disordered configuration `init_disordered.p`, we run `python init_keratin.py` to scale the system and initialise the forces deriving from keratin model. This generates `init_keratin.p` and we then run `python init_steady.py` until it generates a satisfactory initial configuration of keratin distribution. We rename the latest output file to `init_steady.p` and run `python init_cut.py` which cuts a hole in the system and sets the boundary line tension at the hole. The last generated configuration `init_hole.p` can then be simulated with `python -m cells.run -i init_hole.p`.

```python
python init_keratin.py              # scale and initialise keratin
python init_steady.py               # live run to generate initial keratin configuration
mv out.p init_steady.p              # initial keratin configuration
python init_cut.py                  # cut hole and set boundary line tension
python -m cells.run -i init_hole.p  # cutting/closing simulation
```

