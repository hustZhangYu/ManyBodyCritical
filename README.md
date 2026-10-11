# ManyBodyCritical
We consider the many body critical region.\
help from chatgpt

## Source layout and environment

`MbcH.py` builds the Hamiltonian, `JumpOperator.py` builds the spin jump
operators, `RelaxationRate.py` contains perturbation-matrix routines, and
`Main.py` contains sampling experiments. `spin.py` is a local module, not a
package to install from PyPI. Both operator modules use `basis.py` for the
particle-number labels, preserving the original basis order while avoiding
the dense number-operator allocation.

Create an isolated Python environment and install `requirements.txt` with
`python -m pip install -r requirements.txt`. Dependency versions from the
original research environment were not recorded. Run the bounded basis
checks from this directory with `python -m unittest discover -s tests`.

## Research entry points and remaining checks

`Main.py` contains expensive sampling and dense diagonalization; executing
it starts the selected experiment. Check sizes, worker counts and output
CSV names before launching it. The basis refactor does not change the
Hamiltonian, jump operators, sampling parameters or relaxation formulas.

The current `MbcR` combines a fixed-particle-number Hamiltonian with
full-Hilbert-space jump operators. Their dimensions differ for a
nontrivial particle sector. Decide the intended treatment of
particle-changing jumps before using that route; this cleanup does not
silently project them and change the model. Also, `ChangeW14` currently
sets `L=10`, so its name alone is not evidence of a 14-site calculation.
Large sampling runs and scientific outputs have not been reproduced by
the bounded tests.
