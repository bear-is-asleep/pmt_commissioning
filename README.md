# SBND PMT Commissioning

The goal of this work is verify the PMT channel mapping implemented in SBND simulations. We will use a set of east-west filtered cosmics to produce scintillation light in the PMT. We want to predict, only using muon information combined with the detector parameters, the light seen by a given PMT.

$LY = N_\gamma \times f(x) \times \epsilon$

## Efficiency $\epsilon$
Quantum efficiency calculates the PMT's efficiency at detecting the light that does reach it's surface. We have both TPB coated and noncoated PMTs, each of which has a different efficiency depending on $\epsilon_{ws}$ (TPB efficiency). We also have a reflective CPA surface which also is coated with TPB, which means we also have to consider the probability that the light is visible (430 nm) or VUV (128 nm). We consider the probability of a photon reflecting off a TPB coated PMT to be negligble, due to the PMT's high transmissivity. 

We will use binomial probability's to define the probability of a photon being visible:


$P_{vis}(1,c,\epsilon_{ws}) = c \epsilon_{ws} (1-\epsilon_{ws})^{c-1}$ and $P_{VUV} = 1-P_{vis}$


This leads us to the quantum efficiencies for both coated and uncoated PMTs:

$\epsilon_{vis} = QE \times P_{vis}$

$\epsilon_{VUV} = QE\times (\epsilon_{ws} \times P_{vis} + (1-\epsilon_{ws}) \times P_{VUV})$


## Number of photons produced $N_\gamma$
Photons are produced when the minimum ionizing muon passes through the LAr in the TPC, creating VUV photons characterized by the electric field present $N_0$ = 2.4 $\times$ 10$^4$ $\gamma$/MeV [1]. We can calculate the number of photons produced by a track of length $dx$ as follows:

$N_\gamma=N_0\frac{dE}{dx}dx$

Where $\frac{dE}{dx}$ is the muon's MIP energy loss. We assume any PMT that is in the opposite TPC as the track in question will not see any photons. We also receive a number of tracks $N_t$ and a total number of hits $N_h$ (for all tracks) in a single readout. Thus, we split each track length $\ell$ into several steps $dx=\frac{\ell^2}{\ell_{tot}}N_h$ where $\ell_{tot}=\sum_{i=0}^{N_t}\ell_i$, or the total length of all the tracks. Most events only see one cosmic muon, but some see multiple events.

## Fraction of photons reaching the detector
First, we want to calculate the number of photons lost due to ionization, then we want to find how many reach the PMT surface given the reflective surfaces.
### Reflections
Consider a single photon shown in the figure below
### Attenuation $R_I(x)$







